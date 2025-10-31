import argparse
import collections.abc
import itertools
import os
import sys
import typing

import pandas as pd
from tqdm.auto import tqdm
import pysam

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Deduplicate reads in a BAM file based on position and/or tag values.",
        epilog=(
            "Memory usage considerations: (1) if --tag is used, then the maximum number of reads in memory is the "
            "maximum number of reads with the same tag value; (2) if --tag is not used and --record_counts is not "
            "used, then only 1 read is kept in memory at a time, although a dictionary of key values (see --by) from "
            "all unique reads is still used; (3) if --tag is not used but --record_counts is used, then all "
            "deduplicated reads are kept in memory before being written out."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="in.bam",
        help="Path to input BAM file. If not provided, read from standard input."
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="out.bam",
        help="Path to output BAM file. If not provided, write to standard output."
    )
    parser.add_argument(
        "--by",
        required=True,
        help=(
            "Key(s) by which to deduplicate reads. Either a position format ('start' or 'end') or tag name (e.g., 'RX' "
            "for UMIs). Keys can be combined using '&' (AND) or '|' (OR) operators, with '&' operators taking "
            "precedence. For example, 'start&end|RX' will deduplicate reads that share the same start and end "
            "positions, or the same RX tag value. Because no assumptions are made about the order of reads with "
            "respect to their positions in the BAM file, this script may consume significant memory if the --tag "
            "option is not used."
        )
    )
    parser.add_argument(
        "--tag",
        help=(
            "Tag name (e.g., 'CB' for cluster barcodes). Reads must share the same tag value (in addition to any "
            "requirement imposed by --by) to be considered duplicates. If this option is provided, then the input BAM "
            "file is assumed to be sorted by this tag value (e.g., samtools sort -t <tag>). If --by is not a position "
            "format, then this tag must be different than the one specified by --by."
        )
    )
    parser.add_argument(
        "--keep-unmapped",
        action="store_true",
        help=(
            "Keep unmapped reads or reads that do not have the specified tag, as specified by --by. "
            "If not set, these reads will be discarded."
        )
    )
    parser.add_argument(
        "--keep-untagged",
        action="store_true",
        help=(
            "Only relevant if --tag is set. Keep reads that do not have the tag specified by --tag. "
            "If not set, these reads will be discarded."
        )
    )
    parser.add_argument(
        "--record_counts",
        action="store_true",
        help=(
            "Record the number of duplicates represented by each read in a dc tag (analogous to samtools markdup "
            "--duplicate-count). Significantly increase memory usage if the --tag option is not used."
        )
    )
    parser.add_argument(
        '--output_histogram',
        help=(
            "Path to output histogram file; compression is inferred from path extension. Histogram file contains 2 "
            "tab-delimited columns: multiplicity and number of reads with that multiplicity; no header is included."
        )
    )
    parser.add_argument("-u", "--uncompressed", action="store_true", help="Uncompressed BAM output.")
    parser.add_argument("-t", "--threads", type=helpers.positive_int, default=1, help="Number of threads to use.")
    parser.add_argument("--progress", action="store_true", help="Display a progress bar")
    return parser.parse_args()


def deduplicate_reads(
    reads: collections.abc.Iterable,
    bam_out: pysam.AlignmentFile,
    by: str,
    keep_unmapped: bool = False,
    record_counts: bool = False,
    record_histogram: bool = False
) -> tuple[int, int, collections.Counter]:
    """
    Deduplicate reads based on the specified criteria.

    Args
    - reads: Iterable of pysam.AlignedSegment objects (reads).
    - bam_out: pysam.AlignmentFile object for writing deduplicated reads.
    - by: Key(s) by which to deduplicate reads. Either a position format ('start' or 'end') or tag name (e.g., 'RX'
        for UMIs). Keys can be combined using '&' (AND) or '|' (OR) operators, with '&' operators taking precedence. For
        example, 'start&end|RX' will deduplicate reads that share the same start and end positions, or the same RX tag
        value.
    - keep_unmapped: Whether to keep unmapped or untagged reads.
    - record_counts: Whether to record the number of duplicates represented by each read in a dc tag.
    - record_histogram: Whether to record multiplicity histogram.

    Returns
    - n_written: Number of unique reads written to the output BAM.
    - n_total: Total number of reads processed.
    - histogram: (if record_histogram is True) mapping from multiplicity to number of unique reads with that
        multiplicity; (if record_histogram is False) empty Counter.

    Notes
    - Setting `record_counts=True` will increase memory usage significantly, as all primary (non-duplicate) reads must
      be stored in memory until all reads have been checked.
    """
    key_to_idx_primary = dict() # map from "by" value to the index of the first read seen with that value

    track_counts = record_counts or record_histogram
    idx_to_count = collections.Counter() # (only used if track_counts is True) map from index of a read to the number of duplicates represented by that read
    idx_to_read = dict() # (only used if record_counts is True) map from index of a read to the read object itself
    histogram = collections.Counter() # (only used if record_histogram is True) map from multiplicity to number of reads with that multiplicity

    n_written = 0
    idx = -1
    for idx, read in enumerate(reads):
        skip_to_next_read = False
        keys = []
        for tag_set in by.split('|'):
            key = []
            for tag_name in tag_set.split('&'):
                if (tag_name in ('start', 'end') and read.reference_id < 0) or (tag_name not in ('start', 'end') and not read.has_tag(tag_name)):
                    # unmapped read: if keep_unmapped, write it out; otherwise, skip
                    if keep_unmapped:
                        if record_histogram:
                            histogram[1] += 1
                        if record_counts:
                            read.set_tag("dc", 1, value_type="i", replace=True)
                        bam_out.write(read)
                        n_written += 1
                    skip_to_next_read = True
                    break
                if tag_name == 'start':
                    key.extend([read.reference_id, read.reference_start, b's'])
                elif tag_name == "end":
                    key.extend([read.reference_id, read.reference_end, b'e'])
                else:
                    key.append(read.get_tag(tag_name))
            if skip_to_next_read:
                break
            keys.append(tuple(key))
        if skip_to_next_read:
            continue

        primary_read = True
        for key in keys:
            idx_primary = key_to_idx_primary.get(key, None)
            if idx_primary is not None:
                # this read is a duplicate
                primary_read = False
                if track_counts:
                    idx_to_count[idx_primary] += 1
                break
        if primary_read:
            # this read is a primary read
            for key in keys:
                key_to_idx_primary[key] = idx
            if track_counts:
                idx_to_count[idx] += 1
            if record_counts:
                idx_to_read[idx] = read
            else:
                bam_out.write(read)
                n_written += 1

    n_total = idx + 1 # total number of reads processed

    if record_counts:
        for idx, read in idx_to_read.items():
            read.set_tag("dc", idx_to_count[idx], value_type="i", replace=True)
            bam_out.write(read)
            n_written += 1
    if record_histogram:
        histogram.update(collections.Counter(idx_to_count.values()))

    return n_written, n_total, histogram


def deduplicate_bam(
    file_in: str | typing.IO,
    file_out: str | typing.IO,
    by: str,
    tag: str | None = None,
    keep_unmapped: bool = False,
    keep_untagged: bool = False,
    record_counts: bool = False,
    record_histogram: bool = False,
    uncompressed: bool = False,
    threads: int = 1,
    progress: bool = False
) -> tuple[int, int, collections.Counter]:
    """
    Deduplicate reads in a BAM file based on position and/or tag values. The first read in each group of duplicates
    is written to the output BAM file, and subsequent duplicates are discarded.

    Args
    - file_in: Path to input BAM file or a file-like object.
    - file_out: Path to output BAM file or a file-like object.
    - by: Value(s) by which to deduplicate reads. Either a position format ('start', 'end', 'start&end', or 'start|end')
        or tag name (e.g., 'RX' for UMIs).
    - tag: Optional additional tag name to use as criteria for deduplication.
    - keep_unmapped: Whether to keep unmapped or untagged reads, as specified by `by`.
    - keep_untagged: Whether to keep reads without the tag specified by `tag`.
    - record_counts: Whether to record the number of duplicates represented by each read in a dc tag.
    - record_histogram: Whether to record multiplicity histogram.
    - uncompressed: Write uncompressed BAM output.
    - threads: Number of threads to use for writing the output BAM. Only applicable if uncompressed is False.
    - progress: Display a progress bar.

    Returns
    - n_written: Number of unique reads written to the output BAM.
    - n_total: Total number of reads processed.
    - histogram: (if record_histogram is True) mapping from multiplicity to number of unique reads with that
        multiplicity; (if record_histogram is False) empty Counter.
    """
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1
    with pysam.AlignmentFile(file_in, mode="rb", check_sq=False) as bam_in, \
         pysam.AlignmentFile(file_out, mode=mode_out, threads=threads, template=bam_in) as bam_out:
        if tag:
            # first group by tag value (e.g., cluster barcode)
            n_total = 0 # total number of reads processed
            n_written = 0 # number of unique reads written to the output BAM
            histogram = collections.Counter() # map from multiplicity to number of unique reads with that multiplicity
            for tag_value, reads in itertools.groupby(
                tqdm(bam_in.fetch(until_eof=True), disable=not progress),
                key=lambda read: read.get_tag(tag) if read.has_tag(tag) else None
            ):
                if tag_value is None:
                    if keep_untagged:
                        for read in reads:
                            n_total += 1
                            if record_counts:
                                read.set_tag("dc", 1, value_type="i", replace=True)
                            if record_histogram:
                                histogram[1] += 1
                            bam_out.write(read)
                            n_written += 1
                else:
                    cluster_written, cluster_total, cluster_histogram = deduplicate_reads(
                        reads,
                        bam_out,
                        by,
                        keep_unmapped=keep_unmapped,
                        record_counts=record_counts,
                        record_histogram=record_histogram
                    )
                    histogram.update(cluster_histogram)
                    n_total += cluster_total
                    n_written += cluster_written
        else:
            n_written, n_total, histogram = deduplicate_reads(
                tqdm(bam_in.fetch(until_eof=True), disable=not progress),
                bam_out,
                by,
                keep_unmapped=keep_unmapped,
                record_counts=record_counts,
                record_histogram=record_histogram
            )

    if record_histogram:
        # sanity check that histogram sums to total reads written
        assert n_written == sum(k * v for k, v in histogram.items())

    return n_written, n_total, histogram


def main():
    args = parse_args()
    n_written, n_total, histogram = deduplicate_bam(
        args.input if args.input else sys.stdin,
        args.output if args.output else sys.stdout,
        args.by,
        tag=args.tag,
        keep_unmapped=args.keep_unmapped,
        keep_untagged=args.keep_untagged,
        record_counts=args.record_counts,
        record_histogram=bool(args.output_histogram),
        uncompressed=args.uncompressed,
        threads=args.threads,
        progress=args.progress
    )
    if args.output_histogram:
        pd.Series(histogram).sort_index().to_csv(
            args.output_histogram,
            sep='\t',
            index=True,
            header=False
        )
    print(f"Total reads processed: {n_total}. Unique reads written: {n_written}", file=sys.stderr)
    print(f"Duplication rate: {(n_total - n_written) / n_total:.2%}", file=sys.stderr)


if __name__ == "__main__":
    main()
