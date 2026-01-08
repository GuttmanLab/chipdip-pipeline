"""
Count the number of reads at each stage in the ChIP-DIP pipeline.
"""

import argparse
import itertools
import json
import os
import string
import sys
from typing import Union

import numpy as np
import pandas as pd
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count reads at each step of the ChIP-DIP pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "pipeline",
        help=(
            "path to pipeline DAG description file (YAML format). Mapping from level to a mapping describing that "
            "output. Except for a 'data' level, the secondary mapping must contain the key 'path' mapping to a list of "
            "strings that form the path pattern for that output."
        )
    )
    parser.add_argument(
        "--path_samples",
        required=True,
        help="path to samples.json file mapping sample name to mapping from read orientation to list of .fastq.gz files"
    )
    parser.add_argument(
        "--wildcards",
        action='append',
        nargs="+",
        help="wildcards"
    )
    parser.add_argument(
        "--wildcard_names",
        nargs='+',
        help="wildcard names corresponding to --wildcards"
    )
    parser.add_argument(
        "--dir_counts",
        required=True,
        help="path to directory containing count files named {level}[.{wildcard1}.{wildcard2}...].{ext}.count"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="print status messages to stderr"
    )
    parser.add_argument(
        "-o", "--out",
        help="path to save aggregated counts as CSV file"
    )
    parser.add_argument(
        "--dump",
        help="path to dump raw counts as CSV file (columns = level, path, count)"
    )
    parser.add_argument(
        "-s", "--sep",
        default="\t",
        help='separator for printing efficiency output'
    )
    parser.add_argument(
        "--counts_per_template",
        action="store_true",
        help=(
            "report counts per template rather than per read (e.g., for paired-end data, count each pair as one "
            "template). Specifically, use the count_per_template value from the pipeline description file to divide "
            "read counts. Levels without this field are assumed to have count_per_template = 1."
        )
    )
    args = parser.parse_args()
    return args


def path_to_count_ext(path: str) -> str:
    """
    Given a path pattern, return the corresponding count file extension pattern without the leading period.
    """
    if path.endswith(('.fastq.gz', '.fq.gz')):
        return 'fastq-gz.count'
    elif path.endswith('.bam'):
        return 'bam.count'
    else:
        raise ValueError(f"Unknown file extension for path: {path}")


class Node:
    """
    A node represents one step in the ChIP-DIP pipeline for a given sample.
    """

    def __init__(
        self,
        sample: str,
        level: str,
        graph: "Graph",
        children: set[Union["Node", str]] | None = None,
        parent: Union["Node", str, None] = None,
        base: Union["Node", str, None] = None,
        n_reads: int | float = np.nan,
        n_reads_detailed: dict | None = None,
        depth: int | None = None,
    ):
        """
        Args
        - sample: sample name
        - level: level of node in pipeline
        - graph: ChIP-DIP pipeline Graph
        - children: set of children nodes, either as Node objects or node names ("{sample}_{level}")
        - parent: parent node, either as Node object or node name ("{sample}_{level}")
        - base: base node, either as Node object or node name ("{sample}_{level}")
        - n_reads: total number of reads
        - n_reads_detailed: mapping wildcard combinations to number of reads
        - depth: depth of the node in the graph
        """
        if children is None:
            children = set()
        self.sample = sample
        self.level = level
        self.name = f"{sample}_{level}"
        self.graph = graph
        self.children = set()
        self.add_children(children)
        self.set_parent(parent)
        self.base = base
        self.n_reads = n_reads
        self.n_reads_detailed = n_reads_detailed
        if n_reads_detailed:
            n_reads_detailed_sum = np.nansum(list(n_reads_detailed.values()))
            if not np.isnan(n_reads):
                assert n_reads == n_reads_detailed_sum
            else:
                self.n_reads = n_reads_detailed_sum
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
        if isinstance(parent, Node):
            assert parent.graph == self.graph
            self.parent.add_children([self])

    def set_base(self, base):
        self.base = base
        if isinstance(base, Node):
            assert base.graph == self.graph

    def get_root(self):
        node = self
        while node.parent is not None:
            node = node.parent
        return node

    def __repr__(self):
        repr_parent = self.parent.name if isinstance(self.parent, Node) else self.parent
        repr_base = self.base.name if isinstance(self.base, Node) else self.base
        repr_children = (
            ";".join(child.name for child in self.children)
            if len(self.children) > 0
            else self.children
        )
        return f"{self.name} ({self.n_reads} reads, parent={repr_parent}, base={repr_base}, depth={self.depth}, children={repr_children})"

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
            frac_base = self.n_reads / self.base.n_reads if isinstance(self.base, Node) and self.base.n_reads > 0 else np.nan
        frac_parent = None
        if self.parent:
            frac_parent = self.n_reads / self.parent.n_reads if isinstance(self.parent, Node) and self.parent.n_reads > 0 else np.nan
        d = {
            "sample": self.sample,
            "depth": self.depth,
            "level": self.level,
            "parent_level": self.parent.level if isinstance(self.parent, Node) else None,
            "base_level": self.base.level if isinstance(self.base, Node) else None,
            "root_level": root.level if self.parent else None,
            "n_reads": self.n_reads,
            "frac_parent": frac_parent,
            "frac_base": frac_base,
            "frac_root": frac_root,
        }
        return d

    def print_efficiency(self, sep="\t", sample=None, toprint=True):
        d = self.to_dict()
        indent = ""
        if self.depth is not None and self.depth > 0:
            indent = "  " * (self.depth - 1) + "- "
        else:
            if sample is None:
                sample = True
        str_percent_root = (
            f"({100*d['frac_root']:.3g}% of {d['root_level']})"
            if d["frac_root"] is not None
            else ""
        )
        str_percent_base = (
            f"({100*d['frac_base']:.3g}% of {d['base_level']})"
            if d["frac_base"] is not None
            else ""
        )
        str_percent_parent = (
            f"({100*d['frac_parent']:.3g}% of {d['parent_level']})"
            if d["frac_parent"] is not None
            else ""
        )
        percents = sep.join((str_percent_parent, str_percent_base, str_percent_root)).rstrip()

        name = self.name if sample else self.level

        str_efficiency = (
            indent
            + name
            + sep
            + f"{int(self.n_reads) if not np.isnan(self.n_reads) else self.n_reads} reads"
            + ((sep + percents) if percents != "" else "")
        )
        if toprint:
            print(str_efficiency)
        return str_efficiency


class Graph:
    """
    A graph represents a ChIP-DIP pipeline.
    """

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
        return {node.level for node in self.nodes}

    def get_node(self, sample, level):
        matching_nodes = [node for node in self.nodes if node.sample == sample and node.level == level]
        if len(matching_nodes) > 1:
            raise ValueError(f"Multiple nodes found with sample {sample} and level {level}.")
        elif len(matching_nodes) == 0:
            return None
        else:
            return matching_nodes[0]

    def get_node_by_name(self, name):
        if name is None:
            return None
        matching_nodes = [node for node in self.nodes if node.name == name]
        if len(matching_nodes) > 1:
            raise ValueError(f"Multiple nodes found with name {name}.")
        elif len(matching_nodes) == 0:
            return None
        else:
            return matching_nodes[0]

    def get_roots(self):
        return [node for node in self.nodes if node.parent is None]

    def resolve_parents(self):
        """
        Resolve parent nodes given as strings to Node objects.

        Raises: ValueError if a parent cannot be resolved.
        """
        for node in self.nodes:
            if isinstance(node.parent, str):
                parent = self.get_node_by_name(node.parent)
                if isinstance(parent, Node):
                    node.set_parent(parent)
                else:
                    raise ValueError(f"Could not resolve parent {node.parent} for node {node.name}.")

    def resolve_bases(self):
        """
        Resolve base nodes given as strings to Node objects.

        Raises: ValueError if a parent cannot be resolved.
        """
        for node in self.nodes:
            if isinstance(node.base, str):
                base = self.get_node_by_name(node.base)
                if isinstance(base, Node):
                    node.set_base(base)
                else:
                    raise ValueError(f"Could not resolve base {node.base} for node {node.name}.")

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

    def merge_samples_from_node(self, G_merged, node):
        node_parent_level = None if node.parent is None else node.parent.level
        node_merge_parent = (
            G_merged.get_node("all", node_parent_level) if node_parent_level else None
        )
        node_base_level = None if node.base is None else node.base.level
        node_merge_base = G_merged.get_node("all", node_base_level) if node_base_level else None
        existing_merged_levels = G_merged.get_node_levels()
        if node.level not in existing_merged_levels:
            node_merge = Node(
                "all",
                node.level,
                G_merged,
                parent=node_merge_parent,
                base=node_merge_base,
                n_reads=0,
            )
        else:
            node_merge = G_merged.get_node("all", node.level)
            assert node_merge.parent is node_merge_parent
            assert node_merge.base is node_merge_base
        node_merge.n_reads += node.n_reads
        for child in node.children:
            node_child_merge = self.merge_samples_from_node(G_merged, child)
            node_merge.add_children([node_child_merge])
        return node_merge

    def merge_samples(self):
        G_merged = Graph()
        roots = self.get_roots()
        # DFS on G, adding counts and nodes as necessary to G_merged
        for root in roots:
            self.merge_samples_from_node(G_merged, root)
        return G_merged


def collect_pipeline_counts(
    pipeline: dict,
    samples: dict[str, dict[str, list[str]]],
    wildcards: dict[str, list[str]],
    DIR_COUNTS: str,
    count_per_template: bool = False,
    verbose: bool = True,
) -> tuple[Graph, dict[tuple[str, str], int]]:
    """
    Parse pipeline read file counts into a Graph.

    Args
    - pipeline: mapping from level to a structured information (a dictionary) describing that output. Except for the
        "data" level, the info dict must contain the key "path" mapping to a list of strings that form the path pattern
        for that output.
    - samples: mapping of sample names to mapping from read orientation to list of .fastq.gz files
    - wildcards: mapping of placeholder names to lists of values
    - DIR_COUNTS: directory containing count files
    - count_per_template: report counts per template rather than per read (e.g., for paired-end data, count each pair
        as one template). Specifically, use the count_per_template value from the pipeline description file to divide
        read counts. Levels without this field are assumed to have count_per_template = 1
    - verbose: print status messages to standard error.

    Returns
    - G: Graph with potentially multiple root nodes, each corresponding to sample
    """
    G = Graph()
    formatter = string.Formatter()
    all_counts = dict() # map from (level, filename) to count

    for level, node_info in pipeline.items():
        if verbose:
            print(f"Reading counts from {level}", file=sys.stderr)
        scale_factor = node_info.get("count_per_template", 1) if count_per_template else 1
        if level == "data":
            for sample, d in samples.items():
                files_R1 = d["R1"]
                counts = dict()
                for file_number in range(len(files_R1)):
                    path_count = os.path.join(DIR_COUNTS, f'data.{sample}.{file_number}.fastq-gz.count')
                    if not os.path.exists(path_count):
                        counts[(file_number,)] = np.nan
                        print(f"Warning: count file {path_count} does not exist", file=sys.stderr)
                    else:
                        with open(path_count) as f:
                            count = int(f.read().strip()) / scale_factor
                            counts[(file_number,)] = count
                            all_counts[(level, path_count)] = count

                # construct a node for this sample and level; the constructor automatically adds the node to the graph
                Node(
                    sample,
                    level,
                    graph=G,
                    parent=f"{sample}_{node_info['parent']}" if node_info.get("parent") else None,
                    base=f"{sample}_{node_info['base']}" if node_info.get("base") else None,
                    n_reads_detailed=counts
                )
        else:
            fields = sorted(set(
                field for _, field, _, _ in formatter.parse(os.path.join(*node_info['path']))
                if field is not None
            ))
            fields_excluding_sample = [field for field in fields if field != "sample"]
            assert all(field in wildcards.keys() for field in fields_excluding_sample), \
                f"Not all fields {fields} for output level {level} are in wildcards (keys: {wildcards.keys()})."
            path_count_template = os.path.join(
                DIR_COUNTS,
                '.'.join((
                    level,
                    '.'.join([f"{{{field}}}" for field in fields]),
                    path_to_count_ext(node_info["path"][-1])
                ))
            )
            wildcards_filtered = dict()
            for field, values in wildcards.items():
                exclude = set(node_info.get("exclude", dict()).get(field, []))
                wildcards_filtered[field] = [v for v in values if v not in exclude]
            for sample in samples.keys():
                counts = dict()
                for field_combination in itertools.product(*[wildcards_filtered[field] for field in fields_excluding_sample]):
                    path_count = path_count_template.format(**dict(zip(fields_excluding_sample, field_combination), sample=sample))
                    if not os.path.exists(path_count):
                        counts[field_combination] = np.nan
                        print(f"Warning: count file {path_count} does not exist", file=sys.stderr)
                    else:
                        with open(path_count) as f:
                            count = int(f.read().strip()) / scale_factor
                            counts[field_combination] = count
                            all_counts[(level, path_count)] = count

                # construct a node for this sample and level; the constructor automatically adds the node to the graph
                Node(
                    sample,
                    level,
                    graph=G,
                    parent=f"{sample}_{node_info['parent']}" if node_info.get("parent") else None,
                    base=f"{sample}_{node_info['base']}" if node_info.get("base") else None,
                    n_reads_detailed=counts
                )
    G.resolve_bases()
    G.resolve_parents()
    return G, all_counts


def main():
    args = parse_args()

    with open(args.pipeline) as f:
        pipeline = yaml.safe_load(f)

    with open(args.path_samples) as f:
        samples = json.load(f)

    assert len(args.wildcards) == len(args.wildcard_names), \
        "--wildcards and --wildcard_names must have the same number of entries."

    wildcards = dict()
    for i in range(len(args.wildcards)):
        wildcards[args.wildcard_names[i]] = args.wildcards[i]

    G, all_counts = collect_pipeline_counts(
        pipeline,
        samples,
        wildcards,
        args.dir_counts,
        count_per_template=args.counts_per_template,
        verbose=args.verbose
    )
    G_merged = G.merge_samples() if len(samples) > 1 else None

    G.print_efficiency(sep=args.sep, toprint=True)
    if G_merged:
        print()
        G_merged.print_efficiency(sample=False, sep=args.sep, toprint=True)
    if args.out:
        df = G.to_df()
        if G_merged:
            df = pd.concat((df, G_merged.to_df()), axis=0, ignore_index=True)
        df.to_csv(args.out)
    if args.dump:
        pd.Series(all_counts, name='count').rename_axis(['level', 'path']).to_csv(args.dump)


if __name__ == "__main__":
    main()
