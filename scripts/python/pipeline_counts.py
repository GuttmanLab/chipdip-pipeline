'''
Count the number of reads at each stage in the ChIP-DIP pipeline.

Dependencies
- samtools
- standard linux tools: ls, xargs, grep, wc, sort, uniq, awk, sh, gzip
- Python libraries: NumPy, Pandas

Notes
- The pipeline order is specified in the global PIPELINE and PIPELINE_CLUSTERS
  variables.
'''

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import numpy as np
import pandas as pd

regex_fastq = re.compile(r'(\.fastq\.gz$)|(\.fq.gz$)|(\.fastq$)|(\.fq$)')
regex_bam = re.compile(r'(\.bam$)|(\.sam$)|(\.cram$)')
regex_uniq_count = re.compile(r'\s*(?P<count>\d+)\s+(?P<pattern>.*)')

PIPELINE = {
    # (splitfq) raw reads
    'splitfq': {
        'parent': None,
        'path': os.path.join('splitfq', '{aliquot}_R1.part_*.fastq.gz')},

    # (adaptor_trimming_pe) adapter trimming
    'trim': {
        'parent': 'splitfq',
        'path': os.path.join('trimmed', '{aliquot}*val_1.fq.gz')},

    # (barcode_id) barcode idenfication - should not discard any reads; barcode gets appended to read name
    'barcode': {
        'parent': 'trim',
        'path': os.path.join('fastqs', '{aliquot}*_R1*.barcoded.fastq.gz')},

    # (split_bpm_dpm) Split reads into DPM, BPM, short, other
    # - short: any barcode position is [NOT_FOUND]
    # - other: incorrect barcode order
    'short': {
        'parent': 'barcode',
        'path': os.path.join('fastqs', '{aliquot}*_R1*.barcoded_short.fastq.gz')},
    'other': {
        'parent': 'barcode',
        'path': os.path.join('fastqs', '{aliquot}*_R1*.barcoded_other.fastq.gz')},
    'bpm': {
        'parent': 'barcode',
        'path': os.path.join('fastqs', '{aliquot}*_R1*.barcoded_bpm.fastq.gz')},
    'dpm': {
        'parent': 'barcode',
        'path': os.path.join('fastqs', '{aliquot}*_R1*.barcoded_dpm.fastq.gz')},

    ## DNA processing ##
    # (cutadapt_dpm) Remove DPM from read1 of DPM reads
    'dpm-trim': {
        'parent': 'dpm',
        'base': 'dpm',
        'path': os.path.join('trimmed', '{aliquot}*_R1*.barcoded_dpm.RDtrim.fastq.gz')},

    # (bowtie2_align) Align to genome
    'bowtie2': {
        'parent': 'dpm-trim',
        'base': 'dpm',
        'path': os.path.join('alignments_parts', '{aliquot}*.DNA.bowtie2.mapq20.bam')},

    # (add_chr) Chromosome relabeling (add "chr") and filtering (removing non-canonical chromosomes)
    'add-chr': {
        'parent': 'bowtie2',
        'base': 'dpm',
        'path': os.path.join('alignments_parts', '{aliquot}*.DNA.chr.bam')},

    # (repeat_mask)
    'mask': {
        'parent': 'add-chr',
        'base': 'dpm',
        'path': os.path.join('alignments_parts', '{aliquot}*.DNA.chr.masked.bam')},

    # (merge_dna)
    'merge-dna': {
        'parent': 'mask',
        'base': 'dpm',
        'path': os.path.join('alignments', '{aliquot}.DNA.merged.bam')},

    # (thresh_and_split)
    'merge-labeled': {
        'parent': 'merge-dna',
        'base': 'dpm',
        'path': os.path.join('alignments', '{aliquot}.DNA.merged.labeled.bam')},

    # (thresh_and_split)
    'splitbam': {
        'parent': 'merge-labeled',
        'base': 'dpm',
        'path': os.path.join('splitbams', '{aliquot}.DNA.merged.labeled_*.bam')},

    ## Oligo processing ##
    # (cutadapt_oligo)
    'bpm-trim': {
        'parent': 'bpm',
        'base': 'bpm',
        'path': os.path.join('trimmed', '{aliquot}*_R1*.barcoded_bpm.RDtrim.fastq.gz')},

    # (fastq_to_bam)
    'bpm-bam': {
        'parent': 'bpm-trim',
        'base': 'bpm',
        'path': os.path.join('alignments_parts', '{aliquot}*.BPM.bam')},

    # # (merge_beads)
    'bpm-merge': {
        'parent': 'bpm-bam',
        'base': 'bpm',
        'path': os.path.join('alignments', '{aliquot}.merged.BPM.bam')}
}

PIPELINE_CLUSTERS = {
    'cluster': {
        'parent': {'dpm': 'mask', 'bpm': 'bpm-bam'},
        'path': os.path.join('clusters_parts', '{aliquot}*.clusters')},

    'cluster-merge': {
        'parent': {'dpm': 'cluster', 'bpm': 'cluster'},
        'path': os.path.join('clusters', '{aliquot}.clusters')}
}

class Node:
    '''
    A node represents one step in the ChIP-DIP pipeline for a given aliquot.
    '''
    def __init__(self, aliquot, level, graph, children=None, parent=None, base=None, n_reads=np.nan, depth=None):
        '''
        Args
        - aliquot: str
        - level: str
        - graph: Graph
        - children: set of Node or None
        - parent: Node or None
        - base: Node or None
        - n_reads: int or np.nan
        - depth: int or None
        '''
        if children is None:
            children = set()
        self.aliquot = aliquot
        self.level = level
        self.name = f'{aliquot}_{level}'
        self.graph = graph
        self.children = set()
        self.add_children(children)
        self.set_parent(parent)
        self.base = base
        self.n_reads = n_reads
        self.depth = depth
        self.graph.add_node(self)

    def add_children(self, children):
        self.children |= set(children)
        for child in children:
            assert child.graph == self.graph
            assert child.parent in (None, self)
            child.parent = self

    def set_parent(self, parent):
        self.parent = parent
        if parent is not None:
            assert parent.graph == self.graph
            self.parent.add_children([self])

    def get_root(self):
        node = self
        while node.parent is not None:
            node = node.parent
        return node

    def __repr__(self):
        repr_parent = self.parent.name if self.parent is not None else 'None'
        repr_base = self.base.name if self.base is not None else 'None'
        repr_children = ';'.join((child.name for child in self.children)) if len(self.children) > 0 else self.children
        return f'{self.name} ({self.n_reads} reads, parent={repr_parent}, base={repr_base}, depth={self.depth}, children={repr_children})'

    def __hash__(self):
        return hash(self.name)

    def to_dict(self):
        if self.depth is None:
            self.graph.compute_depths()

        frac_root = None
        root = self.get_root()
        if self.parent:
            frac_root = self.n_reads / root.n_reads if root.n_reads > 0 else np.nan

        frac_base = None
        if self.base:
            frac_base = self.n_reads / self.base.n_reads if self.base.n_reads > 0 else np.nan

        frac_parent = None
        if self.parent:
            frac_parent = self.n_reads / self.parent.n_reads if self.parent.n_reads > 0 else np.nan

        d = {
            'aliquot': self.aliquot,
            'depth': self.depth,
            'level': self.level,
            'parent_level': self.parent.level if self.parent is not None else None,
            'base_level': self.base.level if self.base is not None else None,
            'root_level': root.level if self.parent else None,
            'n_reads': self.n_reads,
            'frac_parent': frac_parent,
            'frac_base': frac_base,
            'frac_root': frac_root}
        return d

    def print_efficiency(self, sep='\t', aliquot=None, toprint=True):
        d = self.to_dict()
        indent = ''
        if self.depth > 0:
            indent = '  ' * (self.depth - 1) + '- '
        else:
            if aliquot is None:
                aliquot = True

        str_percent_root = f"({100*d['frac_root']:.3g}% of {d['root_level']})" if d['frac_root'] is not None else ''
        str_percent_base = f"({100*d['frac_base']:.3g}% of {d['base_level']})" if d['frac_base'] is not None else ''
        str_percent_parent = f"({100*d['frac_parent']:.3g}% of {d['parent_level']})" if d['frac_parent'] is not None else ''
        percents = sep.join((str_percent_parent, str_percent_base, str_percent_root)).rstrip()

        name = self.name if aliquot else self.level

        str_efficiency = indent + name + sep + f'{self.n_reads} reads' + ((sep + percents) if percents != '' else '')
        if toprint:
            print(str_efficiency)
        return str_efficiency

class Graph:
    '''
    A graph represents a ChIP-DIP pipeline.
    '''
    def __init__(self, nodes=None):
        if nodes is None:
            nodes = set()
        self.nodes = nodes
        for node in self.nodes:
            assert node.graph in (None, self)
            node.graph = self

    def add_node(self, node):
        if node not in self.nodes:
            assert node.name not in [n.name for n in self.nodes]
            self.nodes.add(node)
            assert node.graph in (None, self)
            node.graph = self

    def get_node_names(self):
        names = [node.name for node in self.nodes]
        names_set = set(names)
        assert len(names_set) == len(names)
        return names_set

    def get_node_levels(self):
        return set([node.level for node in self.nodes])

    def get_node(self, aliquot, level):
        return [node for node in self.nodes if node.aliquot == aliquot and node.level == level][0]

    def get_node_by_name(self, name):
        if name is None:
            return None
        return [node for node in self.nodes if node.name == name][0]

    def get_roots(self):
        return [node for node in self.nodes if node.parent is None]

    def compute_depths_from_node(self, node):
        assert node.depth is not None
        assert node in self.nodes
        for child in node.children:
            assert child.depth in (None, node.depth + 1)
            child.depth = node.depth + 1
            self.compute_depths_from_node(child)

    def compute_depths(self):
        roots = self.get_roots()
        for root in roots:
            assert root.depth in (0, None)
            root.depth = 0
            self.compute_depths_from_node(root)

    def to_df(self, node=None, dicts=None):
        if dicts is None:
            dicts = []
        if node is None:
            next_nodes = self.get_roots()
        else:
            assert node in self.nodes
            dicts.append(node.to_dict())
            next_nodes = node.children
        for next_node in next_nodes:
            self.to_df(node=next_node, dicts=dicts)
        return pd.DataFrame(dicts)

    def print_efficiency(self, node=None, history=None, **kwargs):
        if history is None:
            history = []
        if node is None:
            next_nodes = self.get_roots()
        else:
            assert node in self.nodes
            history.append(node.print_efficiency(**kwargs))
            next_nodes = node.children
        for next_node in next_nodes:
            self.print_efficiency(node=next_node, history=history)
        return history

    def merge_aliquots_from_node(self, G_merged, node):
        node_parent_level = None if node.parent is None else node.parent.level
        node_merge_parent = G_merged.get_node('merged', node_parent_level) if node_parent_level else None
        node_base_level = None if node.base is None else node.base.level
        node_merge_base = G_merged.get_node('merged', node_base_level) if node_base_level else None
        existing_merged_levels = G_merged.get_node_levels()
        if node.level not in existing_merged_levels:
            node_merge = Node('merged', node.level, G_merged, parent=node_merge_parent, base=node_merge_base, n_reads=0)
        else:
            node_merge = G_merged.get_node('merged', node.level)
            assert node_merge.parent is node_merge_parent
            assert node_merge.base is node_merge_base

        node_merge.n_reads += node.n_reads
        for child in node.children:
            node_child_merge = self.merge_aliquots_from_node(G_merged, child)
            node_merge.add_children([node_child_merge])
        return node_merge

    def merge_aliquots(self):
        G_merged = Graph()
        roots = self.get_roots()
        # DFS on G, adding counts and nodes as necessary to G_merged
        for root in roots:
            self.merge_aliquots_from_node(G_merged, root)
        return G_merged

def count_reads(paths, n_processes, unzip=None, filetype=None, agg=True):
    '''
    Count the number of reads in FASTQ or SAM/BAM/CRAM files.

    Args
    - paths: list-like
        Paths to read files. All read files must be of the same file type.
    - n_processes: int
        Number of parallel processes to use.
    - unzip: str. default=None
        Command to write file contents to standard output. If None, defaults to 'gunzip -c'
        - Use 'cat' for uncompressed FASTQ files (i.e., ending in .fastq or .fq)
        - Use 'gunzip -c' for gzip-compressed FASTQ files (i.e., ending in .fastq.gz or .fq.gz)
        - Use 'unpigz -c' for a small speedup for gzip-compressed FASTQ files, if pigz is installed.
    - filetype: str. default=None
        Type of files to count.
        - 'fastq': files ending in .fastq, .fq., .fastq.gz, or .fq.gz
        - 'bam': files ending in .sam, .bam., or .cram
        - None: automatically detect based on the first file in paths
    - agg: bool. default=True
        Return sum of read counts

    Returns: int or list of int
    - If agg is False: returns list of int, giving the read counts for each file in paths
    - If agg is True: returns int, the sum of read counts over all files in paths
    '''
    if unzip is None:
        unzip = 'gunzip -c'
    if filetype is None:
        if regex_fastq.search(paths[0]):
            assert all(map(regex_fastq.search, paths)), 'Files with differing extensions detected.'
            filetype = 'fastq'
        elif regex_bam.search(paths[0]):
            assert all(map(regex_bam.search, paths)), 'Files with differing extensions detected.'
            filetype = 'bam'
        else:
            raise ValueError(f'Unsupported file extension for file {paths[0]}. Must be .fastq.gz, .fq.gz, .bam, .sam., or .cram')
    assert str(filetype).lower() in ('fastq', 'bam')

    paths_concat = ' '.join((f"'{path}'" for path in paths))
    if filetype == 'fastq':
        cmd = f'ls {paths_concat} | xargs -n 1 -P {n_processes} -I % sh -c "{unzip} % | wc -l | awk \'{{print \$1 / 4}}\'"'
    else:
        cmd = f'ls {paths_concat} | xargs -n 1 -P {n_processes} samtools view -@ 4 -c'
    counts = list(map(int, subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip().split()))
    return sum(counts) if agg else counts

def count_reads_clusters(paths, n_processes):
    '''
    Count the number of DPM and BPM reads in SPRITE cluster files.

    Args
    - paths: list-like
        Paths to cluster files
    - n_processes: int
        Number of parallel processes to use.

    Returns: dict
        Keys = 'DPM', 'BPM'
        Values = sum of the number of DPM or BPM reads in cluster files.
    '''
    paths_concat = ' '.join((f"'{path}'" for path in paths))
    cmd = f'LC_ALL=C; ls {paths_concat} | xargs grep -h -F -o -e "DPM" -e "BPM" | sort --parallel={n_processes} | uniq -c'
    counts = dict(DPM=0, BPM=0)
    for line in subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip().split('\n'):
        count, pattern = regex_uniq_count.match(line.strip()).groups()
        counts[pattern] += int(count)
    return counts

def collect_pipeline_counts(
    aliquots,
    DIR_WORKUP,
    n_processes,
    pipeline=None,
    pipeline_clusters=None,
    unzip=None,
    verbose=True):
    '''
    Parse pipeline read file counts into a Graph.

    Args
    - aliquots: list-like of str
    - DIR_WORKUP: str
        Directory to which all paths given in pipeline or pipeline_clusters are relative to
    - n_processes: int
        Number of processes
    - pipeline: dict. default=None
        Pipeline specification. If None, defaults to PIPELINE.
    - pipeline_clusters: dict. default=None
        Pipeline specification for cluster files. If None, defaults to PIPELINE_CLUSTERS.
    - unzip: str. default=None
        Command (e.g., 'zcat', 'unpigz -c') to write file contents to standard output.
        See count_reads() docstring.
    - verbose: bool. default=True
        Print status messages to standard error.

    Returns: tuple(Graph, Graph or None)
    - G: Graph
        Graph with potentially multiple root nodes, each corresponding to an aliquot
    - G_merged: Graph or None
        - If the number of aliquots is 1, G_merged is None
        - Otherwise, a Graph where all the aliquots in G have been merged; one root node
    '''
    if pipeline is None:
        pipeline = PIPELINE
    if pipeline_clusters is None:
        pipeline_clusters = PIPELINE_CLUSTERS
    G = Graph()

    for level, node_info in pipeline.items():
        if verbose:
            print(f'Counting {level}', file=sys.stderr)
        for aliquot in aliquots:
            paths = glob.glob(os.path.join(DIR_WORKUP, node_info['path'].format(aliquot=aliquot)))
            if len(paths) == 0:
                count = 0
            else:
                count = count_reads(paths, n_processes, unzip=unzip)
            base = None
            if 'base' in node_info and node_info['base']:
                base = G.get_node(aliquot, node_info['base'])
            parent = None
            if node_info['parent']:
                parent = G.get_node(aliquot, node_info['parent'])
            node = Node(aliquot, level, graph=G, n_reads=count, parent=parent, base=base)

    for level, node_info in pipeline_clusters.items():
        if verbose:
            print(f'Counting {level}', file=sys.stderr)
        for aliquot in aliquots:
            paths = glob.glob(os.path.join(DIR_WORKUP, node_info['path'].format(aliquot=aliquot)))
            counts = count_reads_clusters(paths, n_processes)
            for base_level, parent in node_info['parent'].items():
                if parent == 'cluster':
                    parent += f'-{base_level}'
                node = Node(
                    aliquot,
                    f'{level}-{base_level}',
                    graph=G,
                    n_reads=counts[base_level.upper()],
                    parent=G.get_node(aliquot, parent),
                    base=G.get_node(aliquot, base_level))

    G_merged = None
    if len(aliquots) > 1:
        G_merged = G.merge_aliquots()
    return G, G_merged

def main():
    parser = argparse.ArgumentParser(description='Count reads at each step of the ChIP-DIP pipeline.')
    parser.add_argument('-d', '--dir-pipeline', help='path to pipeline')
    parser.add_argument('-n', '--n-processes', help='number of processes to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='print status messages to stderr')
    parser.add_argument('-o', '--out', help='path to save counts as CSV file')
    parser.add_argument(
        '-u', '--unzip',
        help='command to write compressed file contents to standard output: e.g., zcat, unpigz -c, etc.')
    parser.add_argument(
        '-s', '--sep',
        default='\t',
        help='separator for printing efficiency output (default="\\t")')

    args = parser.parse_args()
    n_processes = args.n_processes
    DIR_PIPELINE = args.dir_pipeline
    if n_processes is None:
        try:
            n_processes = len(os.sched_getaffinity(0))
        except:
            n_processes = os.cpu_count()

    if DIR_PIPELINE is None:
        DIR_PIPELINE = os.path.abspath(os.getcwd())
    DIR_WORKUP = os.path.join(DIR_PIPELINE, 'workup')

    path_samples_json = os.path.join(DIR_PIPELINE, 'samples.json')
    with open(path_samples_json, 'rt') as f:
        aliquots = tuple(json.load(f).keys())

    G, G_merged = collect_pipeline_counts(
        aliquots,
        DIR_WORKUP,
        n_processes,
        unzip=args.unzip,
        verbose=args.verbose)

    if args.verbose:
        print('Pipline Read Counts Per Aliquot:', file=sys.stderr)
    G.print_efficiency(sep=args.sep)
    if G_merged:
        print()
        if args.verbose:
            print('Pipline Read Counts Aggregated Over All Aliquots', file=sys.stderr)
        G_merged.print_efficiency(aliquot=False, sep=args.sep)
    if args.out:
        df = G.to_df()
        if G_merged:
            df = pd.concat((df, G_merged.to_df()), axis=0, ignore_index=True)
        df.to_csv(args.out)

if __name__ == '__main__':
    main()
