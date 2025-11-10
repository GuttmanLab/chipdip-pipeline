import argparse
import contextlib
import itertools
import os
import sys
import typing

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers

import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collapse or decollapse SAM/BAM files with paired reads. When collapsing paired reads, the output is a 27 "
            "(or 28 if --chrom_map is provided) column BED-like format with the following fields: <chrom> <start> "
            "<end> <11 required fields of the first read> <TLV-encoded tags of the first read> <11 required fields of "
            "the second read> <TLV-encoded tags of the second read>. Both reads in a pair must be aligned to the same "
            "chromosome; otherwise, both reads are discarded. Tags are encoded in a type-length-value format so that "
            "the tab character can be safely used as a field separator in the BED-like format. TLV encoding is "
            "performed as follows: [length]:[tag]:[type]:[value] [length]:[tag]:[type]:[value] ... where 'length' is "
            "the number of characters of the tag field, and each tag field (previous tab-delimited) is now delimited "
            "by a single space character. For example, the tag 'NM:i:1' would be encoded as '6:NM:i:1', since the "
            "length of the string 'NM:i:1' is 6 characters."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        help="Path to name-collated SAM/BAM file containing paired reads. If not provided, read from standard input.",
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to output file. If not provided, write to standard output.",
    )
    parser.add_argument(
        "--decollapse",
        action="store_true",
        help="Decollapse a BED-like file back into a SAM/BAM file with paired reads.",
    )
    parser.add_argument(
        "--keep-unpaired",
        action="store_true",
        help=(
            "Retain unpaired reads in the output. In the collapsed BED-like format, the non-existant 'pair' of the "
            "unpaired read will have its read name set to '*'."
        ),
    )
    parser.add_argument(
        "--chrom_map",
        help=(
            "Chromosome name map file. If provided, prepend an extra field with the 0-based index of the chromosome "
            "to each output line (e.g., for subsequent sorting by a custom chromosome order)."
        )
    )
    parser.add_argument(
        '--template_bam',
        help="(decollapse mode) Path to a template BAM file from which to copy the header for the output BAM."
    )
    parser.add_argument(
        "-u", "--uncompressed",
        action="store_true",
        help="(If --decollpase is used) Output uncompressed BAM."
    )
    parser.add_argument(
        "--threads",
        type=helpers.positive_int,
        default=1,
        help="(If --decollapse is used without --uncompressed) Number of threads to use for BAM compression."
    )
    return parser.parse_args()


@contextlib.contextmanager
def default_file_object(path=None, mode='wt', default=None):
    '''
    If path is not None, open it with the provided mode within a with statement.
    If path is None, yield default. (If default is None, yields None.)

    Args
    - path: None | str
        Path to file
    - mode: str. default='wt'
        Mode with which to open path.
    - default: None | file object. default=None.
        Default file object to yield if path is None.
        Typically, set to sys.stdin, sys.stdout or sys.stderr.

    Yields: file object or None

    Example
        with default_file_object(path, default=sys.stderr) as f:
            f.write('hi')
    '''
    if path is None:
        yield default
    else:
        with open(path, mode=mode) as file:
            yield file


def encode_tags(tags) -> str:
    """Encode SAM optional fields as concatenated <len>:<field> TLV blocks."""
    out = []
    for field in tags:
        out.append(f"{len(field)}:{field}")
    return " ".join(out)


def decode_tags(tlv_str: str) -> list[str]:
    """
    Decode TLV-encoded SAM optional fields into a list of strings.

    Returns: list of fields (e.g., ['NM:i:1', 'MD:Z:35A4', ...])
    - Empty list if tlv_str is empty.
    """
    tags = []
    i = 0
    n = len(tlv_str)
    while i < n:
        # get length of field
        j = tlv_str.find(":", i)
        if j == -1:
            raise ValueError(f"Error decoding TLV string: could not find ':' after position {i} in '{tlv_str}'.")
        length_str = tlv_str[i:j]
        try:
            length = int(length_str)
        except ValueError:
            raise ValueError(f"Error decoding TLV string: could not parse length '{length_str}' as integer.")

        # Move to start of field; skip ':'
        j += 1
        end_field = j + length
        if end_field > n:
            raise ValueError(f"Error decoding TLV string: expected field of length {length} at position {j}, "
                             f"but only {n - j} characters remain in '{tlv_str}'.")
        field = tlv_str[j:end_field]
        tags.append(field)

        # Move to next block, skip space
        i = end_field + 1
    return tags


def collapse_sam(
    file_in: str | typing.IO,
    file_out: str | None = None,
    path_chrom_map: str | None = None,
    keep_unpaired: bool = False,
):
    """Collapse paired reads in a SAM/BAM file into a BED-like format."""
    chrom_to_index = None
    if path_chrom_map is not None:
        chrom_map = helpers.parse_chrom_map(path_chrom_map)
        chrom_to_index = {v: str(i) for i, v in enumerate(chrom_map.values())}
    with pysam.AlignmentFile(file_in, mode="r") as bam_in, \
         default_file_object(file_out, default=sys.stdout) as f_out:
        for name, pair in itertools.groupby(
            bam_in.fetch(until_eof=True),
            key=lambda read: read.query_name
        ):
            pair = list(pair)
            if keep_unpaired or len(pair) == 2:
                if len(pair) == 1:
                    mate1 = pair[0]
                    mate2 = None
                    start = str(mate1.reference_start)
                    end = str(mate1.reference_end)
                else:
                    mate1, mate2 = pair
                    assert mate1.reference_id == mate2.reference_id, (
                        f"Error: Read pair '{name}' aligned to different chromosomes: {mate1.reference_name} and "
                        f"{mate2.reference_name}"
                    )
                    start = str(min(mate1.reference_start, mate2.reference_start))
                    end = str(max(mate1.reference_end, mate2.reference_end))
                chrom = str(mate1.reference_name)
                fields = mate1.to_string().split('\t')
                tag_str = encode_tags(fields[11:]) if len(fields) > 11 else ''
                mate1_str = '\t'.join(fields[:11] + [tag_str])
                if mate2 is not None:
                    fields = mate2.to_string().split('\t')
                    tag_str = encode_tags(fields[11:]) if len(fields) > 11 else ''
                    mate2_str = '\t'.join(fields[:11] + [tag_str])
                else:
                    mate2_str = '*' + '\t' * 11
                if chrom_to_index:
                    out_str = '\t'.join([chrom_to_index[chrom], chrom, start, end]) + '\t' + mate1_str + '\t' + mate2_str
                else:
                    out_str = '\t'.join([chrom, start, end]) + '\t' + mate1_str + '\t' + mate2_str
                f_out.write(out_str + '\n')


def decollapse_bam(
    file_in: str | typing.IO,
    file_out: str | typing.IO,
    path_template_bam: str,
    keep_unpaired: bool = False,
    uncompressed: bool = False,
    threads: int = 1,
):
    """
    Decollapse a BED-like file back into a SAM/BAM file with paired reads.

    Args
    - file_in: path or file object to input BED-like file
        Columns: 1-3 = chrom, start, end
                 4-14 = 11 required fields of first read
                 15 = TLV-encoded tags of first read
                 16-26 = 11 required fields of second read
                 27 = TLV-encoded tags of second read
    - file_out: path or file object to output BAM file
    """
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1
    with pysam.AlignmentFile(path_template_bam, "r") as f:
        header_dict = f.header.to_dict()
    with default_file_object(file_in, mode='rt', default=sys.stdin) as f_in, \
         pysam.AlignmentFile(file_out, mode=mode_out, header=header_dict, threads=threads) as bam_out:
        header = bam_out.header
        for line in f_in:
            if not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            assert len(cols) == 27, f"Error: Expected 27 columns in collapsed BED-like format, found {len(cols)}."
            if keep_unpaired or cols[15] != '*':
                mate1_tags = decode_tags(cols[14])
                read = pysam.AlignedSegment.fromstring('\t'.join(cols[3:14] + mate1_tags), header)
                bam_out.write(read)
            if cols[15] != '*': # mate 2 exists and is paired
                mate2_tags = decode_tags(cols[26])
                read = pysam.AlignedSegment.fromstring('\t'.join(cols[15:26] + mate2_tags), header)
                bam_out.write(read)


def main():
    args = parse_args()

    if args.decollapse:
        decollapse_bam(
            args.input,
            args.output,
            path_template_bam=args.template_bam,
            keep_unpaired=args.keep_unpaired,
            uncompressed=args.uncompressed,
            threads=args.threads,
        )
    else:
        collapse_sam(
            args.input if args.input else sys.stdin,
            args.output,
            keep_unpaired=args.keep_unpaired,
            path_chrom_map=args.chrom_map,
        )


if __name__ == "__main__":
    main()