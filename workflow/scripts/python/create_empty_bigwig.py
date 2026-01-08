import argparse

import pyBigWig
import pysam


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create an empty bigWig file.")
    parser.add_argument("output", help="Path to output empty bigWig file")

    # require either BAM file or chromosome sizes file to generate bigWig file header
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--bam', help="Path to SAM/BAM file to extract chromosome sizes")
    group.add_argument('--chroms', help="Path to chromosome sizes file (tab-separated: chrom, size)")

    return parser.parse_args()


def get_chrom_sizes_from_bam(path: str) -> dict[str, int]:
    """Extract chromosome sizes from SAM/BAM file header"""
    with pysam.AlignmentFile(path, "r") as f:
        header = f.header
        chrom_sizes = dict(zip(header.references, header.lengths))
    return chrom_sizes


def parse_chrom_sizes(path) -> dict[str, int]:
    """Parse chromosome sizes file (tab-separated: chrom, size)"""
    chrom_sizes = dict()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                chrom, size = line.split('\t')
                chrom_sizes[chrom] = int(size)
    return chrom_sizes


def write_empty_bigwig(path, chrom_sizes: dict[str, int]) -> None:
    with pyBigWig.open(path, "w") as bw:
        bw.addHeader(list(chrom_sizes.items()))


def main():
    args = parse_arguments()

    if args.bam:
        chrom_sizes = get_chrom_sizes_from_bam(args.bam)
    else:
        chrom_sizes = parse_chrom_sizes(args.chroms)

    write_empty_bigwig(args.output, chrom_sizes)


if __name__ == "__main__":
    main()
