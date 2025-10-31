import argparse
import os
import re
import sys
import typing

import pysam

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert an antibody oligo FASTQ file to a BAM file."
    )
    parser.add_argument(
        "input",
        metavar="in.fastq",
        help=(
            "Path to input FASTQ file with barcodes in the read name (e.g., read::[tag1][tag2][tag3]) "
            "and read sequences starting with a UMI."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="out.bam",
        help=(
            "Path to output BAM file. If not provided, write to stdout. Reads will be unaligned but have SAM/BAM tags "
            "for read type (RT:Z:<str>), barcode (CB:Z:<str>), and UMI (RX:Z:<sequence>; optionally QX:Z:<qualities>)."
        )
    )
    parser.add_argument("--UMI_len", metavar="#", type=int, required=True, help="Length of the UMI")
    parser.add_argument("--num_tags",
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
        "--remove_UMI_from_seq",
        action="store_true",
        help="Remove UMI from the sequence and quality fields after it is copied to the tags.",
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
    parser.add_argument("--drop_quals", action="store_true", help="Do not store FASTQ quality scores in the BAM file.")
    parser.add_argument(
        "--read_type_prefixes",
        nargs="+",
        default=["BEAD_"],
        help="Prefixes to identify read types in the barcode. Reads with these prefixes will be processed.",
    )
    parser.add_argument("-u", "--uncompressed", action="store_true", help="Uncompressed BAM output.")
    parser.add_argument("-t", "--threads", type=helpers.positive_int, default=1, help="Number of threads to use.")
    parser.add_argument("--verbose", action="store_true", help="Print progress to stderr.")
    return parser.parse_args()


def convert_reads(
    path_in: str,
    path_out: str | typing.IO,
    num_tags: int,
    UMI_length: int,
    header: pysam.AlignmentHeader | dict | None = None,
    add_sample_to_barcode: str | None = None,
    read_type_prefixes: tuple[str, ...] = ("BEAD_",),
    remove_barcode_from_names: bool = False,
    remove_UMI_from_seq: bool = False,
    tag_read_type: str = "RT",
    tag_barcode: str = "CB",
    drop_quals: bool = False,
    uncompressed: bool = False,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, int]:
    """
    Convert an antibody oligo FASTQ file to a BAM file.

    Args
    - path_in: Path to input FASTQ file.
    - path_out: Path to output BAM file, or a file object to write to.
    - num_tags: Number of tags in barcode,
    - UMI_length: Length of the UMI.
    - header: SAM/BAM format header. If None, defaults to "@HD	VN:1.6	SO:unknown".
    - add_sample_to_barcode: Sample/aliquot name to add to the barcode.
    - read_type_prefixes: Prefixes to identify read types in the barcode.
        Reads with these prefixes will be processed. Otherwise, they will be skipped.
    - remove_barcode_from_names: Remove barcode from read names.
    - remove_UMI_from_seq: Remove UMI from read sequences and qualities.
    - tag_read_type: 2-letter SAM/BAM tag name for the read type.
    - tag_barcode: 2-letter SAM/BAM tag name for the cluster barcode.
    - drop_quals: Do not store FASTQ quality scores in the BAM file.
    - uncompressed: Write uncompressed BAM output.
    - threads: Number of threads to use for writing the output BAM. Only applicable if uncompressed is False.
    - verbose: Print progress to stderr.

    Returns: dictionary with the following keys:
    - "input": Number of reads processed.
    - "written": Number of read written.
    - "skipped": Number of reads skipped due to errors (i.e., no barcode detected in read name).
    """
    if header is None:
        header = dict(HD=dict(VN='1.6', SO='unknown'))
    barcode_suffix = f'.{add_sample_to_barcode}' if add_sample_to_barcode is not None else ''
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1

    # regular expression for bead oligo name
    PATTERN = re.compile("::" + num_tags * r"\[([a-zA-Z0-9_-]+)\]")

    n_total, n_written, n_skipped = 0, 0, 0
    with pysam.AlignmentFile(path_out, mode=mode_out, header=header, threads=threads) as f_out, \
         helpers.file_open(path_in, mode="rt") as f_in:
        for qname, seq, _, quals in helpers.fastq_parse(f_in):
            n_total += 1
            if verbose and n_total % 100000 == 0:
                print('Reads processed:', n_total, file=sys.stderr)
            match = PATTERN.search(qname)
            if match:
                read_type = match.group(1)
                if read_type.startswith(read_type_prefixes):
                    barcode = '.'.join(match.groups()[1:]) + barcode_suffix

                    # create read object
                    aligned_segment = pysam.AlignedSegment()
                    aligned_segment.query_name = (PATTERN.sub("", qname) if remove_barcode_from_names else qname).lstrip('@')
                    aligned_segment.query_sequence = seq[UMI_length:] if remove_UMI_from_seq else seq
                    aligned_segment.flag = 4 # unmapped
                    if not drop_quals:
                        aligned_segment.query_qualities_str = quals[UMI_length:] if remove_UMI_from_seq else quals
                        aligned_segment.set_tag("QX", quals[:UMI_length], value_type="Z", replace=True)
                    aligned_segment.set_tag(tag_read_type, read_type, value_type="Z", replace=True)
                    aligned_segment.set_tag(tag_barcode, barcode, value_type="Z", replace=True)
                    aligned_segment.set_tag("RX", seq[:UMI_length], value_type="Z", replace=True)

                    f_out.write(aligned_segment)
                    n_written += 1
                else:
                    print(
                        f"Warning: No read type of {read_type_prefixes} detected in first barcode position of",
                        f"read '{qname}'. Skipping read.",
                        file=sys.stderr
                    )
                    n_skipped += 1
            else:
                n_skipped += 1
                if verbose:
                    print(
                        f"Warning: No barcode detected in read name '{qname}'. Skipping read.",
                        file=sys.stderr
                    )
            if verbose and n_total % 100000 == 0:
                print('Reads processed:', n_total, file=sys.stderr)

    reads_dict = dict(input=n_total, written=n_written, skipped=n_skipped)
    if verbose:
        print('Read counts:', reads_dict, file=sys.stderr)
    return reads_dict


def main() -> None:
    args = parse_arguments()
    convert_reads(
        args.input,
        args.output if args.output else sys.stdout,
        args.num_tags,
        args.UMI_len,
        add_sample_to_barcode=args.add_sample_to_barcode,
        read_type_prefixes=tuple(args.read_type_prefixes),
        remove_barcode_from_names=args.remove_barcode_from_names,
        remove_UMI_from_seq=args.remove_UMI_from_seq,
        tag_read_type=args.tag_read_type,
        tag_barcode=args.tag_barcode,
        drop_quals=args.drop_quals,
        uncompressed=args.uncompressed,
        threads=args.threads,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
