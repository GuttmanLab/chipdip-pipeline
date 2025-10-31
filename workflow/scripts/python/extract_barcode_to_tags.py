import argparse
import os
import re
import typing
import sys

import pysam

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers

'''
Extract barcode from read name to SAM/BAM read type and barcode tags.

Example input read:
read1::[DPM1][Y10][E4][O62][E31][O47]	16	chr1	10015	30	95M	*	0	0	*	*

Example output read (default options):
read1::[DPM1][Y10][E4][O62][E31][O47]	16	chr1	10015	30	95M	*	0	0	*	*	RT:Z:DPM1	CB:Z:Y10.E4.O62.E31.O47

Example output read (options --remove_barcode_from_names --add_sample_to_barcode aliquot1 --reverse_barcode_order):
read1	16	chr1	10015	30	95M	*	0	0	*	*	RT:Z:DPM1	CB:Z:O47.E31.O62.E4.Y10.aliquot1
'''

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract barcode from read name to SAM/BAM read type and barcode tags."
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="in.bam",
        help=(
            "Path to input BAM file containing reads with barcodes in read names. "
            "If not provided, read from standard input."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="out.bam",
        help="Path to output BAM file with read type and barcode tags. If not provided, write to standard output.",
    )
    parser.add_argument(
        "--num_tags",
        type=int,
        required=True,
        help="Number of tags in barcode, including DPM/BEAD tags.",
    )
    parser.add_argument(
        "--add_sample_to_barcode",
        metavar="sample_name",
        help="Add a sample/aliquot name to the barcode."
    )
    parser.add_argument(
        "-r",
        "--remove_barcode_from_names",
        action="store_true",
        help="Remove the barcode from read names after it is copied to the tags.",
    )
    parser.add_argument(
        '--reverse_barcode_order',
        action='store_true',
        help="Reverse the order of barcode tags when adding to the barcode tag.",
    )
    parser.add_argument(
        '--add_tags',
        type=str,
        action='append',
        help="Add custom SAM/BAM tags to all reads. TAG:TYPE:VALUE."
    )
    parser.add_argument(
        "--tag_read_type",
        default="RT",
        help="2-letter SAM/BAM tag name for the read type.",
    )
    parser.add_argument(
        "--tag_barcode",
        default="CB",
        help="2-letter SAM/BAM tag name for the cluster barcode. Barcodes will be formatted as tag1.tag2[.sample]",
    )
    parser.add_argument(
        "--read_type_prefixes",
        nargs="+",
        default=["DPM", "BEAD_"],
        help="Prefixes to identify read types in the barcode. Reads with these prefixes will be processed.",
    )
    parser.add_argument("-u", "--uncompressed", action="store_true", help="Uncompressed BAM output.")
    parser.add_argument("-t", "--threads", type=helpers.positive_int, default=1, help="Number of threads to use.")
    parser.add_argument("--verbose", action="store_true", help="Print progress to stderr.")
    return parser.parse_args()


def extract_barcode(
    file_in: str | typing.IO,
    file_out: str | typing.IO,
    num_tags: int,
    add_sample_to_barcode: str | None = None,
    read_type_prefixes: tuple[str, ...] = ("DPM", "BEAD_"),
    remove_barcode_from_names: bool = True,
    reverse_barcode_order: bool = False,
    add_tags: list[str] | None = None,
    tag_read_type: str = "RT",
    tag_barcode: str = "CB",
    uncompressed: bool = False,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, int]:
    """
    Extract barcode from read name to SAM/BAM read type and barcode tags.

    Args
    - file_in: Path to BAM file containing reads with barcodes in read names.
    - file_out: Path to output BAM file with read type and barcode tags.
    - num_tags: Number of tags in barcode.
    - add_sample_to_barcode: Sample/aliquot name to add to the barcode.
    - read_type_prefixes: Prefixes to identify read types in the barcode.
        Reads with these prefixes will be processed. Otherwise, they will be skipped.
    - remove_barcode_from_names: Remove barcode from read names.
    - reverse_barcode_order: Reverse the order of barcode tags when adding to the barcode tag.
    - add_tags: Custom SAM/BAM tags (TAG:TYPE:VALUE) to add to all reads.
    - tag_read_type: 2-letter SAM/BAM tag name for the read type.
    - tag_barcode: 2-letter SAM/BAM tag name for the cluster barcode.
    - uncompressed: Write uncompressed BAM output.
    - threads: Number of threads to use for writing the output BAM. Only applicable if uncompressed is False.
    - verbose: Print progress to stderr.

    Returns: dictionary with the following keys:
    - "input": Number of reads processed
    - "written": Number of read written
    - "skipped": Number of reads skipped due to errors (i.e., no barcode detected in read name)
    """
    barcode_suffix = f'.{add_sample_to_barcode}' if add_sample_to_barcode is not None else ''
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1
    add_tags_parts = [] # list of [TAG, TYPE, VALUE]
    if add_tags is not None:
        for tag in add_tags:
            tag_parts = tag.split(':', maxsplit=2)
            if len(tag_parts) != 3:
                raise ValueError(f"Error: Invalid tag format '{tag}'. Expected TAG:TYPE:VALUE.")
            if len(tag_parts[0]) != 2:
                raise ValueError(f"Error: Invalid tag name '{tag_parts[0]}'. Expected 2-letter tag name.")
            add_tags_parts.append(tag_parts)
    pattern = re.compile("::" + num_tags * r"\[([a-zA-Z0-9_-]+)\]")
    n_total, n_written, n_skipped = 0, 0, 0
    with pysam.AlignmentFile(file_in, mode="r") as bam_in, \
         pysam.AlignmentFile(file_out, mode=mode_out, threads=threads, template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            n_total += 1
            qname_copy = read.query_name
            match = pattern.search(qname_copy)
            if match:
                if remove_barcode_from_names:
                    read.query_name = pattern.sub("", qname_copy)
                read_type = match.group(1)
                if read_type.startswith(read_type_prefixes):
                    if reverse_barcode_order:
                        barcode = '.'.join(match.groups()[:0:-1]) + barcode_suffix
                    else:
                        barcode = '.'.join(match.groups()[1:]) + barcode_suffix
                    read.set_tag(tag_read_type, read_type, value_type="Z", replace=True)
                    read.set_tag(tag_barcode, barcode, value_type="Z", replace=True)
                    for tag_name, tag_type, value in add_tags_parts:
                        read.set_tag(tag_name, value, value_type=tag_type, replace=True)
                    bam_out.write(read)
                    n_written += 1
                else:
                    print(
                        f"Warning: No read type of {read_type_prefixes} detected in first barcode position of",
                        f"read '{qname_copy}'. Skipping read.",
                        file=sys.stderr
                    )
                    n_skipped += 1
            else:
                n_skipped += 1
                if verbose:
                    print(
                        f"Warning: No barcode detected in read name '{qname_copy}'. Skipping read.",
                        file=sys.stderr
                    )
            if verbose and n_total % 100000 == 0:
                print('Reads processed:', n_total, file=sys.stderr)

    reads_dict = dict(input=n_total, written=n_written, skipped=n_skipped)
    if verbose:
        print('Read counts:', reads_dict, file=sys.stderr)
    return reads_dict


def main():
    args = parse_args()
    extract_barcode(
        args.input if args.input else sys.stdin,
        args.output if args.output else sys.stdout,
        args.num_tags,
        add_sample_to_barcode=args.add_sample_to_barcode,
        read_type_prefixes=tuple(args.read_type_prefixes),
        remove_barcode_from_names=args.remove_barcode_from_names,
        reverse_barcode_order=args.reverse_barcode_order,
        add_tags=args.add_tags,
        tag_read_type=args.tag_read_type,
        tag_barcode=args.tag_barcode,
        uncompressed=args.uncompressed,
        threads=args.threads,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
