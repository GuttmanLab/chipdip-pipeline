"""
Calculate effective genome size from BAM file header and excluding masked regions
"""

import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import file_open, positive_int, parse_chrom_map
import pysam


def main():
    args = parse_arguments()
    size = effective_genome_size(args.input, args.mask, args.chrom_map, args.threads)
    print(size)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate effective genome size from BAM file header and excluding masked regions."
    )
    parser.add_argument("input", metavar="in.bam", help="Input BAM file")
    parser.add_argument("-m", "--mask", metavar="PATH", help="Mask file")
    parser.add_argument("-c", "--chrom_map", metavar="PATH", help="Chromosome name map file")
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        default=1,
        metavar="#",
        help="Number of threads to use for compressing/decompressing BAM files",
    )
    return parser.parse_args()


def effective_genome_size(path_bam, path_mask=None, path_chrom_map=None, threads=1):
    """
    Effective genome size = size of genome - size of mask, filtered for selected chromosomes

    Args
    - path_bam: str
        Path to BAM file
    - path_mask: str. default=None
        Path to a merged, coordinate-sorted BED file
        - e.g., sort -k1,1 -k2,2n in.bed | bedtools merge
        Supports uncompressed and gzip-compressed files.
        If None, mask size is assumed to be 0.
    - path_chrom_map: str. default=None
        Path to a chromosome name map file.
        If None, all chromosomes in the BAM file header are allowed.
    - threads: int. default=1
        Number of threads to use for compressing/decompressing BAM files

    Returns
    - effective_genome_size: int
        Number of base pairs
    """
    # select set of chromosomes: intersection of chromosome names in BAM file header and chrom_map
    with pysam.AlignmentFile(path_bam, "rb", threads=threads) as f:
        ref_seqs = f.header.to_dict()["SQ"]
    chrom_set = set(ref_seq["SN"] for ref_seq in ref_seqs)
    if path_chrom_map is not None:
        chrom_map = parse_chrom_map(path_chrom_map)
        chrom_set &= set(chrom_map.keys())

    # calculate genome size of BAM file
    genome_size = sum(
        ref_seq["LN"] for ref_seq in ref_seqs if chrom_set is None or ref_seq["SN"] in chrom_set
    )

    # calculate mask size
    if path_mask is not None:
        mask_size = calculate_bed_size(path_mask, chrom_set)
    else:
        mask_size = 0

    return genome_size - mask_size


def calculate_bed_size(path_bed, chrom_set=None):
    """
    Calculate length spanned by regions in a BED file, excluding specified chromosomes.

    Args
    - path_mask: str
        Path to a merged, coordinate-sorted BED file
        - e.g., sort -k1,1 -k2,2n in.bed | bedtools merge
        Supports uncompressed and gzip-compressed files.
    - chrom_set: set of str. default=None
        If provided, only include regions from these chromosomes.

    Returns
    - size: int
        Number of base pairs
    """
    size = 0
    with file_open(path_bed) as f:
        for line in f:
            chrom, start, end = line.decode().strip().split()[:3]
            if chrom_set is None or chrom in chrom_set:
                size += int(end) - int(start)
    return size


if __name__ == "__main__":
    main()
