import cluster as c
import argparse
import pysam
import re
import os
from collections import defaultdict


def parse_args():

    parser = argparse.ArgumentParser(description='Add protein label to DNA bamfile')
    parser.add_argument('-i', '--input_bam', 
                        dest='input_bam', 
                        type=str, 
                        required=True,
                        help='Aligned DNA Bamfile')
    parser.add_argument('-o', '--output_bam', 
                        dest='output_bam', 
                        type=str, 
                        required=True,
                        help='Path to output bam with protein tag added')
    parser.add_argument('-c', '--clusters', 
                        dest='clusters', 
                        type=str, 
                        required=True, 
                        help='Clusters from which to assign protein tag')
    parser.add_argument('--num_tags', 
                        dest='num_tags', 
                        type=int, 
                        required=True, 
                        help='Number of tags in barcode')
    parser.add_argument('-l', '--labeled', 
                        dest='labeled', 
                        type=bool, 
                        action='store',
                        required=True, 
                        help='Does cluster file have clusters labeled?')
    parser.add_argument('--min_oligos',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        required=False,
                        help = "The minimum number of oligos to call a cluster")
    parser.add_argument('-t', '--threshold',
                        action = "store",
                        type = float,
                        required = False,
                        help = "The fraction of oligos on a bead needed to label the bead.")
    parser.add_argument('--max_size',
                        action = "store",
                        type = int,
                        required = False,
                        help = "The maximum cluster size to keep")

    opts = parser.parse_args()

    return opts

def main():
    args = parse_args()
    if args.labeled==True:
        labels = load_labels(args.clusters)
    else:
        labels = call_labels(args.clusters,  args.threshold, args.min_oligos, args.max_size)
    label_bam_file(args.input_bam, args.output_bam, labels, args.num_tags)


def load_labels(clusterfile):
    labels = {}
    with open(clusterfile, 'r') as clusters:
        for line in clusters:
            labeled_barcode = line.split('\t', 1)[0]
            label, barcode = labeled_barcode.split('.', 1)
            labels[barcode] = label
    return labels

def call_labels(clusterfile, threshold, min_oligos, max_size):
    labels = {}
    with open(clusterfile, 'r') as clusters:
        for line in clusters:
            barcode, *reads = line.rstrip('\n').split('\t')
            label = c.label_cluster_reads(reads, threshold, min_oligos, max_size)
            labels[barcode] = label
    return labels

def label_bam_file(input_bam, output_bam, labels, num_tags):
    count, skipped, duplicates, written = 0, 0, 0, 0
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    library = os.path.basename(input_bam).split('.')[0]
    header = construct_read_group_header(input_bam, labels)
    found = defaultdict(set)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for read in in_bam.fetch(until_eof = True):
            count += 1
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())[1:]
            barcode.append(library)
            full_barcode = '.'.join(barcode)
            position = read.reference_name + "_" + str(read.reference_start)
            if position in found[full_barcode]:
                duplicates +=1
            else:
                label = labels[full_barcode]
                found[full_barcode].add(position)
                try:
                    read.tags += [('RG', label)]
                    out_bam.write(read)
                    written +=1
                except KeyError:
                    skipped += 1

    print('Total reads:', count)
    print('Reads written:', written)
    print('Duplicate reads:', duplicates)
    print('Reads with an error not written out:', skipped)

def construct_read_group_header(input_bam, labels):
    proteins = set(labels.values())
    with pysam.AlignmentFile(input_bam, 'rb') as input_file:
        bam_header = input_file.header.to_dict()
        read_group_dict = [{"ID":name} for name in list(proteins)]
        bam_header["RG"] = read_group_dict
        return bam_header

if __name__== "__main__":
    main()
