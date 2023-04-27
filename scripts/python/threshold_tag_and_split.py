import argparse
import pysam
import re
import os
from collections import defaultdict, Counter
from pathlib import Path
import tqdm

"""
Generate individual bam files for each antibody from a master bam file of all read alignments based on cluster assignments. Reads are deduplicated as individual files are generated (ie. only one read per start position per cluster is kept).
"""


def parse_args():

    parser = argparse.ArgumentParser(
        description="Add antibody label to DNA bamfile and generate individual bamfiles for each antibody"
    )
    parser.add_argument(
        "-i",
        "--input_bam",
        dest="input_bam",
        type=str,
        required=True,
        help="Master aligned DNA Bamfile",
    )
    parser.add_argument(
        "-o",
        "--output_bam",
        dest="output_bam",
        type=str,
        required=True,
        help="Path to output master bam with antibody tag added",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        dest="clusters",
        type=str,
        required=True,
        help="Path to clusters used to assign antibody labels",
    )
    parser.add_argument(
        "--num_tags",
        dest="num_tags",
        type=int,
        required=True,
        help="Number of tags in barcode",
    )
    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        type=str,
        action="store",
        required=True,
        help="Directory to write individual antibody bams to",
    )
    parser.add_argument(
        "--min_oligos",
        action="store",
        type=int,
        required=True,
        help="The minimum number of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--proportion",
        action="store",
        type=float,
        required=True,
        help="The maximum representation proportion of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--max_size",
        action="store",
        type=int,
        required=True,
        help="The maximum cluster size to keep",
    )

    opts = parser.parse_args()

    return opts


def main():
    args = parse_args()
    print("Using min_oligos: ", args.min_oligos)
    print("Using proportion: ", args.proportion)
    print("Using max_size: ", args.max_size)
    print("Writing bam to: ", args.output_bam)
    print("Writing splitbams to: ", args.dir)
    labels = assign_labels(
        args.clusters, args.min_oligos, args.proportion, args.max_size
    )
    label_bam_file(args.input_bam, args.output_bam, labels, args.num_tags)
    split_bam_by_RG(args.output_bam, args.dir)


def assign_labels(clusterfile, min_oligos, proportion, max_size):
    """
    Assign labels to all clusters based on oligo reads

    Args:
        clusterfile(str): Path to clusterfile
        min_oligos(int): number of oligos that are of one type that need to be exceeded to assign the cluster (ie. stricktly greater than)
        proportion(float): fraction of oligo reads that are of one type needed to assign the cluster
        max_size(int): maximum dna reads allowed per cluster
    """
    labels = {}
    with open(clusterfile, "r") as clusters:
        for line in tqdm.tqdm(clusters):
            barcode, *reads = line.rstrip("\n").split("\t")
            multilabel = label_cluster_reads(reads, min_oligos, proportion, max_size)
            labels[barcode] = multilabel
    return labels


def label_cluster_reads(reads, min_oligos, threshold, max_size):
    """
    Assign a label to a cluster based on oligo reads

    Args:
        reads(list): list of cluster formated reads
        threshold(float): fraction of oligo reads that are of one type needed to assign the cluster
        min_oligos(int): number of oligos of one type that need to be exceeded to assign the cluster (ie. strictly greater than)
        max_size(int): maximum dna reads allowed per cluster
    """
    bead_reads = [read for read in reads if read.startswith("BPM")]
    if len(bead_reads) == 0:
        return "none"
    cluster_size = len(reads) - len(bead_reads)
    if int(cluster_size) > int(max_size):
        return "filtered"
    bead_labels = Counter([read.split(":")[0].split("_", 1)[1] for read in bead_reads])
    candidate = bead_labels.most_common()[0]
    if candidate[1] < min_oligos:
        return "uncertain"
    elif candidate[1] / sum(bead_labels.values()) < threshold:
        return "ambiguous"
    else:
        return candidate[0]
    return "malformed"


def label_bam_file(input_bam, output_bam, labels, num_tags):
    """
    Add antibody label to individual reads of the master DNA bam file

    Args:
        input_bam(str): Path to input master bam file
        output_bam(str): Path to write labeled bam file
        labels(dict): Dictionary of cluster assignments [barcode-> antibody]
        num_tags(int): number of tags in barcode
    """
    count, duplicates, skipped = 0, 0, 0
    written = defaultdict(int)
    pattern = re.compile("::" + num_tags * "\[([a-zA-Z0-9_\-]+)\]")
    library = os.path.basename(input_bam).split(".")[0]
    header = construct_read_group_header(input_bam, labels)
    found = defaultdict(set)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for read in in_bam.fetch(until_eof=True):
            count += 1
            if count % 10000000 == 0:
                print(count)
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())[1:]
            barcode.append(library)
            full_barcode = ".".join(barcode)
            position = read.reference_name + "_" + str(read.reference_start)
            if position in found[full_barcode]:
                duplicates += 1
            else:
                try:
                    readlabel = labels[full_barcode]
                    found[full_barcode].add(position)
                    read.set_tag("RG", readlabel, replace=True)
                    out_bam.write(read)
                    written[readlabel] += 1
                except KeyError:
                    skipped += 1

    print("Total reads:", count)
    print("Reads written:", written)
    print("Duplicate reads:", duplicates)
    print("Reads with an error not written out:", skipped)


def construct_read_group_header(input_bam, labels):
    """
    Add read group tags for each antibody to the bam file header

    Args:
        input_bam(str): Path to input master bam file
        labels(dict): Dictionary of cluster assignments [barcode-> antibody]
    """
    proteins = set(labels.values())
    sample_name = input_bam.split(".", 1)[0]
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        bam_header = input_file.header.to_dict()
        read_group_dict = [{"ID": name, "SM": sample_name} for name in list(proteins)]
        bam_header["RG"] = read_group_dict
        return bam_header


def split_bam_by_RG(output_bam, output_dir):
    """
    Split a bam file based on the read group

    Args:
        output_bam(str): Path to bam file containing reads with read group assignments
        output_dir(str): Directory to write read group specific bam files to
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    format_string = output_dir + "/%*_%!.%."
    pysam.split("-f", format_string, output_bam)


if __name__ == "__main__":
    main()
