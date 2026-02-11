import argparse
import collections
import contextlib
import typing
import sys

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.figure
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Remove unpaired reads."
    )
    parser.add_argument(
        "-i", "--input",
        metavar="in.bam",
        help=(
            "Input BAM file, with reads grouped by name (i.e., after running samtools collate) such that reads from a "
            "read pair are adjacent. If not provided, read from standard input."
        )
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-o", "--output",
        metavar="out.bam",
        help="Output BAM file for paired primary alignments. If neither --output nor --no-output is set, write to standard out."
    )
    group.add_argument(
        "--no_output",
        action="store_true",
        help=(
            "Do not write paired reads to output. Still allow other outputs (unpaired reads and fragment length "
            "distribution)."
        )
    )
    parser.add_argument(
        "--unpaired",
        metavar="unpaired.bam",
        help="Optional output BAM file for unpaired, unmapped, secondary, or supplementary reads."
    )
    parser.add_argument(
        "-l", "--length_dist",
        metavar="length_dist.tsv",
        help="Output 2-column (fragment_length, count) tab-delimited histogram file of fragment length distribution."
    )
    parser.add_argument(
        "-p", "--length_dist_plot",
        metavar="length_dist_plot.pdf",
        help="Output histogram plot of fragment length distribution."
    )
    parser.add_argument(
        "--recalc_tlen",
        action="store_true",
        help=(
            "Re-calculate template length instead of relying on the TLEN field in the BAM file. This re-calculated "
            "template length value follows the current SAM specification (v1.6) for TLEN, which is the length from the "
            "leftmost mapped base to the rightmost mapped base of a read pair, excluding soft-clipped bases. "
            "Example use case: When bowtie2 (as of version 2.5.4) is run in local alignment mode without the "
            "--soft-clipped-unmapped-tlen flag, the TLEN value includes soft-clipped bases. "
            "Note: This calculation assumes that primary read pairs are in FR orientation (reads on opposite strands "
            "and pointing towards each other; see https://www.htslib.org/algorithms/duplicate.html). Filtering the "
            "input BAM file with samtools view -f 3 (require that each read is paired (1) and mate mapped in proper "
            "pair (2)) is typically sufficient. However, since the 0x2 FLAG is aligner-dependent, one could use "
            "samtools view -e '((flag & 0x30) == 0x10 || (flag & 0x30) == 0x20) && ((flag & 0xC0) == 0x40 || (flag & "
            "0xC0) == 0x80) && flag.paired' to ensure that the read is marked as read1 or read2, but not both, and "
            "either the read or its mate is mapped to the reverse strand, but not both."
        )
    )
    parser.add_argument(
        "-u", "--uncompressed",
        action="store_true",
        help="(Only relevant if --no_output is not set) Uncompressed paired reads BAM output."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        metavar="#",
        help="Number of threads for compressing/decompressing BAM files",
    )
    return parser.parse_args()


def remove_unpaired(
    file_in: str | typing.IO,
    file_out: str | typing.IO | None = None,
    file_unpaired: str | typing.IO | None = None,
    recalc_tlen: bool = False,
    uncompressed: bool = False,
    threads: int = 1,
) -> tuple[collections.Counter[int], int]:
    '''
    Remove unpaired reads by read name. Unmapped, secondary, or supplementary reads are also considered "unpaired."

    Args
    - file_in: path to name-collated BAM file. If None, read from standard in.
    - file_out: path to output BAM file of paired reads. If None, no output is written (use if only want histogram).
    - file_unpaired: optional path to output BAM file of unpaired reads. If None, do not write unpaired reads.
    - recalc_tlen: whether to re-calculate template length instead of relying on the TLEN field in the BAM file.
        This re-calculated template length value follows the current SAM specification (v1.6) for TLEN, which is the
        length from the leftmost mapped base to the rightmost mapped base of a read pair, excluding soft-clipped bases.
        Example use case: When bowtie2 (as of version 2.5.4) is run in local alignment mode without the
        --soft-clipped-unmapped-tlen flag, the TLEN value includes soft-clipped bases.
        Note: This calculation assumes that primary read pairs are in FR orientation (reads on opposite strands and
        pointing towards each other; see https://www.htslib.org/algorithms/duplicate.html). See the command-line
        argument description for --recalc_tlen for example an samtools view filtering command to ensure that this
        assumption is met.
    - uncompressed: Whether to write uncompressed BAM output.
    - threads: Number of threads to use for reading and writing BAM files

    Returns
    - c: Counter mapping from fragment lengths to number of paired reads
    - n_unpaired: number of unpaired reads
    '''
    in_pair = False
    paired_read = None
    c = collections.Counter()
    n_unpaired = 0

    with contextlib.ExitStack() as stack:
        bam_in = stack.enter_context(pysam.AlignmentFile(file_in, 'rb', threads=threads))
        header = bam_in.header.to_dict()
        if file_out:
            bam_out = stack.enter_context(pysam.AlignmentFile(
                file_out,
                mode="wbu" if uncompressed else "wb",
                threads=1 if uncompressed else threads,
                header=header
            ))
        else:
            bam_out = None
        if file_unpaired:
            bam_unpaired = stack.enter_context(pysam.AlignmentFile(file_unpaired, 'wb', threads=threads, header=header))
        else:
            bam_unpaired = None

        for read in bam_in.fetch(until_eof=True):
            # Skip unpaired, unmapped, secondary, or supplementary reads
            if (not read.is_paired) or read.is_unmapped or read.is_secondary or read.is_supplementary:
                n_unpaired += 1
                if bam_unpaired:
                    bam_unpaired.write(read)
                continue

            if in_pair and read.qname == paired_read.qname:
                # check that alignment is not chimeric; removing supplementary reads above should already ensure this
                assert paired_read.reference_id == read.reference_id

                if recalc_tlen:
                    # check that read pair is in FR orientation
                    assert (read.is_forward and paired_read.is_reverse) or (paired_read.is_forward and read.is_reverse)

                    # Calculate template length based on current SAM specification (v1.6; commit 94500cf), which is the
                    # length from the leftmost mapped base to the rightmost mapped base of a read pair, excluding
                    # soft-clipped bases. This means that for dove-tailed read pairs, template length is calculated as
                    # "TLEN#2" as shown in footnote 16 of the SAM specification.
                    start = min(paired_read.reference_start, read.reference_start)
                    end = max(paired_read.reference_end, read.reference_end)
                    template_length = end - start
                    assert template_length > 0
                    if template_length != abs(read.template_length):
                        print(
                            (f"Warning: Re-calculated template length {template_length} does not match TLEN "
                             f"{read.template_length} for read {read.qname}."),
                            file=sys.stderr
                        )
                else:
                    assert paired_read.template_length == -read.template_length
                    template_length = abs(read.template_length)
                c[template_length] += 1
                if bam_out:
                    bam_out.write(paired_read)
                    bam_out.write(read)
                in_pair = False
            else:
                if in_pair:
                    n_unpaired += 1
                    if bam_unpaired:
                        bam_unpaired.write(paired_read)
                paired_read = read
                in_pair = True
        if in_pair:
            n_unpaired += 1
            if bam_unpaired:
                bam_unpaired.write(paired_read)
    return c, n_unpaired


def plot_length_distribution(c: collections.Counter[int], path_out: str) -> matplotlib.figure.Figure:
    '''
    Plot histogram of fragment length distribution.

    Args
    - c: Counter of fragment lengths
    - path_out: Path to output histogram plot file

    Returns: matplotlib figure
    '''
    lengths = list(c.keys())
    counts = list(c.values())

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.bar(lengths, counts, width=1.0, color='blue')
    ax.set_xlabel('Fragment Length')
    ax.set_ylabel('Number of Paired Reads')
    ax.set_title('Fragment Length Distribution')
    if len(c) > 0:
        ax.set_xticks(range(0, max(lengths) + 1, 100))
    else:
        ax.set_xticks(range(0, 1000, 1))
    fig.savefig(path_out)
    return fig


def main():
    args = parse_arguments()

    if args.output:
        file_out = args.output
    elif args.no_output:
        file_out = None
        if not (args.unpaired or args.length_dist or args.length_dist_plot):
            print("Error: no outputs are specified.", file=sys.stderr)
            sys.exit(1)
    else:
        file_out = sys.stdout

    c, n_unpaired = remove_unpaired(
        args.input if args.input else sys.stdin,
        file_out=file_out,
        file_unpaired=args.unpaired,
        recalc_tlen=args.recalc_tlen,
        uncompressed=args.uncompressed,
        threads=args.threads
    )

    if args.length_dist:
        df = pd.DataFrame({'fragment_length': list(c.keys()), 'count': list(c.values())})
        df.sort_values('fragment_length').to_csv(args.length_dist, sep='\t', index=False)
    if args.length_dist_plot:
        plot_length_distribution(c, args.length_dist_plot)
        plt.close()

    print(f"Number of reads that are not primary paired reads: {n_unpaired}", file=sys.stderr)


if __name__ == '__main__':
    main()