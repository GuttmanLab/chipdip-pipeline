import argparse
import collections
import contextlib
import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import pandas as pd
import matplotlib.pyplot as plt
import pysam


def main():
    args = parse_arguments()
    c, n_unpaired = remove_unpaired(
        args.input,
        path_out_bam=args.output,
        path_length_dist=args.length_dist,
        local=args.local,
        threads=args.threads
    )
    print(f"Number of unpaired reads: {n_unpaired}", file=sys.stderr)
    if args.length_dist_plot:
        plot_length_distribution(c, args.length_dist_plot)


def remove_unpaired(
    path_in_bam: str | None = None,
    path_out_bam: str | None = None,
    path_unpaired_bam: str | None = None,
    path_length_dist: str | None = None,
    local: bool = False,
    threads: int = 1,
) -> tuple[collections.Counter[int], int]:
    '''
    Remove unpaired reads by read name.

    Args
    - path_in_bam: path to name-collated BAM file. If None, read from standard in.
    - path_out_bam: path to output BAM file of paired reads. If None, write to standard out.
    - path_unpaired_bam: optional path to output BAM file of unpaired reads. If None, do not write unpaired reads.
    - path_length_dist: optional path to histogram of fragment length distribution as a 2-column tab-delimited table,
        where the first column is fragment length, and the second column is the number of paired reads with that length.
    - local: wWhether the input BAM file contains local alignments. If True, do not include soft-clipped bases in 
        fragment length. Manually calculate fragment lengths instead of using the TLEN SAM field.
    - threads: Number of threads to use for reading and writing BAM files

    Returns
    - c: Counter mapping from fragment lengths to number of paired reads
    - n_unpaired: number of unpaired reads
    '''
    path_out_bam = path_out_bam if path_out_bam is not None else sys.stdout
    path_in_bam = path_in_bam if path_in_bam is not None else sys.stdin

    in_pair = False
    paired_read = None
    c = collections.Counter()
    n_unpaired = 0
    with contextlib.ExitStack() as stack:
        file_in = stack.enter_context(pysam.AlignmentFile(path_in_bam, 'rb', threads=threads))
    # with pysam.AlignmentFile(path_in_bam, 'rb', threads=threads) as file_in:
        header = file_in.header.to_dict()
        file_out = stack.enter_context(pysam.AlignmentFile(path_out_bam, 'wb', threads=threads, header=header))
        # with pysam.AlignmentFile(path_out_bam, 'wb', threads=threads, header=header) as file_out:
        if path_unpaired_bam:
            file_unpaired = stack.enter_context(pysam.AlignmentFile(path_unpaired_bam, 'wb', threads=threads, header=header))
        else:
            file_unpaired = None
        for read in file_in.fetch(until_eof=True):
            if in_pair and read.qname == paired_read.qname:
                assert paired_read.reference_name == read.reference_name
                assert paired_read.reference_end >= paired_read.reference_start
                assert read.reference_end >= read.reference_start
                assert paired_read.reference_id == read.reference_id
                assert paired_read.template_length == -read.template_length
                if paired_read.is_reverse:
                    assert read.is_forward
                    template_length = paired_read.reference_end - read.reference_start
                    if not local:
                        assert template_length == read.template_length
                else:
                    assert paired_read.is_forward and read.is_reverse
                    template_length = read.reference_end - paired_read.reference_start
                    if not local:
                        assert template_length == paired_read.template_length
                c[template_length] += 1
                file_out.write(paired_read)
                file_out.write(read)
                in_pair = False
            else:
                if in_pair:
                    n_unpaired += 1
                    if file_unpaired:
                        file_unpaired.write(paired_read)
                paired_read = read
                in_pair = True
    if path_length_dist:
        s = pd.Series(c, name='num_paired_reads')
        s.index.name = 'template_length'
        s.sort_index().to_csv(
            path_length_dist,
            sep='\t',
            header=True,
            index=True
        )
    return c, n_unpaired


def plot_length_distribution(c: collections.Counter[int], path_out: str) -> None:
    '''
    Plot histogram of fragment length distribution.

    Args
    - c: Counter of fragment lengths
    - path_out: Path to output histogram plot file
    '''
    lengths = list(c.keys())
    counts = list(c.values())

    plt.figure(figsize=(10, 6))
    plt.bar(lengths, counts, width=1.0, color='blue')
    plt.xlabel('Fragment Length')
    plt.ylabel('Number of Paired Reads')
    plt.title('Fragment Length Distribution')
    if len(c) > 0:
        plt.xticks(range(0, max(lengths) + 1, 100))
    else:
        plt.xticks(range(0, 1000, 1))
    plt.savefig(path_out)
    plt.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Remove unpaired reads based on identical genomic alignment coordinates."
    )
    parser.add_argument(
        "-i", "--input",
        metavar="in.bam",
        help=("Input BAM file, with reads grouped by name (i.e., after running "
              "samtools collate) such that reads from a read pair are adjacent. "
              "If not provided, read from standard input.")
    )
    parser.add_argument(
        "-o", "--output",
        metavar="out.bam",
        help="Output BAM file. If not provided, write to standard out."
    )
    parser.add_argument(
        "-l", "--length-dist",
        metavar="length_dist.tsv",
        help="Output 2-column tab-delimited file for histogram of fragment length distribution. "
    )
    parser.add_argument(
        "-p", "--length-dist-plot",
        metavar="length_dist_plot.pdf",
        help="Output histogram plot of fragment length distribution."
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help=(
            "Whether the input BAM file contains local alignments. If True, do not include soft-clipped bases in "
            "fragment length. Manually calculate fragment lengths instead of using the TLEN SAM field."
        )
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        metavar="#",
        help="Number of threads for compressing/decompressing BAM files",
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()