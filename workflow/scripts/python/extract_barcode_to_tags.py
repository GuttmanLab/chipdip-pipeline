import argparse
import collections.abc
import contextlib
import json
import os
import re
import sys
import typing

from tqdm.auto import tqdm
import pysam

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers

REGEX_READ_NAME = re.compile(r'^[>@]?(?P<name>\S+?)::(?P<barcode>(?:\[[^\]]+\])+)(?:\s+(?P<sam_tags>.+))?$')

SAM_TAG_ARRAY_TO_PYTHON_TYPE = dict(
    c=int, # int8
    C=int, # uint8
    s=int, # int16
    S=int, # uint16
    i=int, # int32
    I=int, # uint32
    f=float # float
) # see SAM specification version 1.6, section 1.5 (https://samtools.github.io/hts-specs/SAMv1.pdf)

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert an antibody oligo FASTA/FASTQ file to a BAM file of unmapped reads.",
        epilog=(
            "Example usage: python bpm_fastx_to_bam.py -i in.fastq -o out.bam --num_tags 1 4 1 "
            "--add_sample_to_barcode aliquot1 --remove_barcode_from_names --read_type_prefixes BEAD_"
            "--discard_inconsistent_R1_R2 --require_read_type. "
            "Example input read name: '@r1::[BEAD_1234][BEAD_1234][TAG1][TAG2][TAG3] RX:Z:AA' "
            "Example output SAM entry: 'r1	4	*	0	255	*	*	0	0	*	*	RT:Z:BEAD_1234	"
            "CB:Z:TAG1.TAG2.TAG3.aliquot1	RX:Z:AA"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="in.fasta|in.fastq",
        help=(
            "Path to input file. Assumed to be FASTQ unless --fasta set. If not provided, read from standard input. "
            "Read names are prefixed with '>' or '@' and are formatted as '<name>::<barcode>[ <SAM_tags>]' where "
            "<name> is a SAM-format-compatible read name, "
            "<barcode> is a series of bracketed tags (e.g., [tag1][tag2][tag3]) with barcode tags detected in the R1 "
            "read preceding barcode tags detected in the R2 read (such as output by splitcode --mod-names), and "
            "<SAM_tags> are optional whitespace-delimited SAM tags (e.g., RX:Z:[UMI1]-[UMI2])."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="out.bam",
        help="Path to output BAM file. If not provided, write to stdout."
    )
    parser.add_argument(
        "-d",
        "--discarded",
        metavar="discarded.bam",
        help="Path to output BAM file for discarded reads. An XD tag is added indicating the reason for discard."
    )
    parser.add_argument(
        "--input_fmt",
        choices=('fasta', 'fastq', 'bam'),
        default='fastq',
        help="Input file format."
    )
    parser.add_argument(
        "--num_tags",
        type=int,
        nargs=3,
        metavar=("NUM_TAGS_IN_R1", "NUM_TAGS_IN_R2", "NUM_TAGS_OVERLAP"),
        required=True,
        help=(
            "Number of barcode tags contributed by each read orientation. NUM_TAGS_OVERLAP specifies the number of "
            "tags sequenced from both orientations and must be <= min(NUM_TAGS_IN_R1, NUM_TAGS_IN_R2). Barcode tags "
            "in the read name are assumed to be ordered 5'->3' from each read orientation. For example, if 2 tags are "
            "sequenced in R1 and 4 tags are sequenced in R2, with 1 overlapping tag, then the barcode would look like "
            "[tag1][tag2][tag5][tag4][tag3][tag2], and --num_tags should be 2 4 1. Overlapping tags are used to check "
            "for consistency. The final barcode is taken from the NUM_TAGS_IN_R1 tags from the R1 read and the first "
            "(NUM_TAGS_IN_R2 - NUM_TAGS_OVERLAP) tags from the R2 read. If NUM_TAGS_TAGS_IN_R1 + NUM_TAGS_IN_R2 - "
            "NUM_TAGS_OVERLAP is less than the total number of tags in the barcode, then the extra tags at the end of "
            "the barcode are ignored."
        )
    )

    # read processing options
    parser.add_argument(
        "--remove_barcode_from_names",
        action="store_true",
        help="Remove the barcode from read names after it is copied to the barcode and read type tags.",
    )
    parser.add_argument(
        '--reverse_barcode_order',
        action='store_true',
        help="Reverse the order of R2 barcode tags when adding to the SAM/BAM barcode tag.",
    )
    parser.add_argument(
        "--add_sample_to_barcode",
        metavar="sample_name",
        help="Add a sample name to the barcode."
    )
    parser.add_argument(
        '--add_SAM_tags',
        metavar='TAG:TYPE:VALUE',
        action='append',
        help=(
            "Add custom SAM/BAM tags to all reads. Takes precedence over SAM tags extracted from the read name or "
            "--extract_UMI."
        )
    )
    parser.add_argument(
        "--tag_read_type",
        default="RT",
        help="2-letter SAM/BAM tag name for the read type.",
    )
    parser.add_argument(
        "--tag_barcode",
        default="CB",
        help="2-letter SAM/BAM tag name for the barcode. Barcodes will be formatted as tag1.tag2[.sample]",
    )
    parser.add_argument(
        "--read_type_prefixes",
        action='append',
        help=(
            "Prefixes to identify read type tags in the barcode. Read type tags are removed from the barcode and "
            "stored in a separate read type SAM tag. Multiple read type tag values are joined with a hyphen."
        )
    )
    parser.add_argument(
        "--extract_UMI",
        metavar="start:stop:#",
        help=(
            "Location of sequence to extract and add to RX SAM tag. start:stop are 0-based (Pythonic) indices. "
            "# must be either 0 (remove UMI from sequence after extraction) or 1 (keep UMI in sequence). "
            "If input is FASTQ format, then corresponding quality scores are added to a QX SAM tag unless --drop_quals "
            "is set. If an RX or QX SAM tag already exists in the read name, then the extracted UMI and quality scores "
            "are appended with a hyphen (e.g., RX:Z:[UMI1]-[UMI2])."
        )
    )

    # read filtering options
    parser.add_argument(
        "--discard_inconsistent_R1_R2",
        action='store_true',
        help=(
            "Check that tags sequenced from both R1 and R2 orientations were identified consistently by comparing tag "
            "names; discard read pairs if any inconsistency."
        )
    )
    parser.add_argument(
        '--require_read_type',
        action='store_true',
        help=(
            "(If --read_type_prefixes is not passed an empty string) Require that one of the read type prefixes be "
            "present in the barcode; otherwise, discard the read."
        )
    )
    parser.add_argument(
        '--discard_UMI_N',
        action='store_true',
        help=(
            "Discard reads with 'N' in the UMI sequence (whether extracted via --extract_UMI or in a pre-existing RX "
            "tag)."
        )
    )
    parser.add_argument(
        '--discard_UMI_mismatch',
        action='store_true',
        help='Discard reads where UMI sequences extracted from different sources do not match. Implies --discard_UMI_N'
    )

    # output options
    parser.add_argument(
        "--drop_quals",
        action="store_true",
        help="(Only for FASTQ input) Do not store FASTQ quality scores in the BAM file."
    )
    parser.add_argument(
        "-u", "--uncompressed",
        action="store_true",
        help="Uncompressed BAM output."
    )
    parser.add_argument(
        "--output_counts",
        metavar="counts.json",
        help=(
            "Output summary counts JSON file. Keys = total, written, r1_r2_mismatch, no_read_type, umi_mismatch. "
            "Keys r1_r2_mismatch, no_read_type, umi_mismatch correspond to the number of reads matching the filters "
            "for --discard_inconsistent_R1_R2, --require_read_type, and (--discard_UMI_N or --discard_UMI_mismatch), "
            "respectively. r1_r2_mismatch is null if NUM_TAGS_OVERLAP is 0. no_read_type is null if "
            "--read_type_prefixes is not specified. For these 3 keys, positive values indicate the number of reads "
            "written out; negative values indicate discarded reads, such that written = total + sum of negative "
            "counts. In other words, any discarded read is only counted in one category."
        )
    )

    # performance options
    parser.add_argument(
        "-t", "--threads",
        type=helpers.positive_int,
        default=1,
        help="Only relevant if -u/--uncompressed is not set. Number of BAM output compression threads to use."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress to stderr."
    )
    return parser.parse_args()


@contextlib.contextmanager
def open_reads_file(file: typing.IO | str, mode: str, fmt: str = 'fastq', **kwargs) -> collections.abc.Generator[typing.IO]:
    '''
    Transparently open a file or yield a file object within a with statement.

    Args
    - file: path to a file, or a file object
        If a file object is given (e.g., sys.stdout), then the file object is NOT closed after use.
    - mode: file object mode
    - fmt: If 'SAM' or 'BAM' (not case sensitive), use pysam.AlignmentFile to open the file, then yield the file within
        a with statement.
    - **kwargs: additional keyword arguments to pass to pysam.AlignmentFile() if fmt is 'SAM' or 'BAM'

    Yields: file object or None
    '''
    if fmt.lower() in ('sam', 'bam'):
        with pysam.AlignmentFile(file, mode=mode, **kwargs) as f:
            yield f
    elif isinstance(file, str):
        with helpers.file_open(file, mode=mode) as f:
            yield f
    else:
        yield file


def convert_SAM_tag_value(
    value_type: typing.Literal["A", "i", "f", "Z", "H", "B"],
    value_str: str
) -> str | int | float | list[int | float]:
    """
    Convert a SAM tag value string to the appropriate Python type based on the value type.

    Args
    - value_type: SAM tag value type
    - value_str: SAM tag value string

    Returns: SAM tag value in appropriate Python type
    - 'Z', 'H', 'A' -> str
    - 'i' -> int
    - 'f' -> float
    - 'B' -> list of int or float
    """
    if value_type == 'A': # printable character
        assert len(value_str) == 1, f"Error: Invalid BAM char tag value '{value_str}'."
        return value_str
    elif value_type == 'i': # signed integer
        return int(value_str)
    elif value_type == 'f': # single-precision float
        return float(value_str)
    elif value_type in ('Z', 'H'):
        # Z can be any printable string, including space, but not tab
        # H is a hex string; I'll let pysam perform any necessary validation when actually adding the tag to the read
        # via pysam.AlignedSegment.set_tag()
        return value_str
    elif value_type == 'B': # integer or numeric array
        value_parts = value_str.split(',')
        if len(value_parts) < 2:
            raise ValueError(f"Error: Invalid BAM array tag value '{value_str}'.")
        array_type = value_parts[0]
        if array_type not in SAM_TAG_ARRAY_TO_PYTHON_TYPE:
            raise ValueError(f"Error: Invalid BAM array tag type '{array_type}'.")
        return [SAM_TAG_ARRAY_TO_PYTHON_TYPE[array_type](x) for x in value_parts[1:]]
    else:
        raise ValueError(f"Error: Invalid SAM tag type '{value_type}'.")


def process_SAM_tag_str(s: str) -> tuple[str, str, str | int | float | list[int | float]]:
    """
    Process a SAM tag string of the form TAG:TYPE:VALUE into (TAG, TYPE, VALUE), where VALUE is converted to the
    appropriate Python type, such as for use with pysam.AlignedSegment.set_tag(TAG, VALUE, value_type=TYPE).
    """
    tag, value_type, value = s.split(':', maxsplit=2)
    if len(tag) != 2:
        raise ValueError(f"Error: Invalid tag name '{tag}'. Expected 2-letter tag name.")
    value = convert_SAM_tag_value(value_type, value)
    return tag, value_type, value


def process_SAM_tags(
    sam_tags_list: list[tuple[str, str, str | int | float | list[int | float]]],
    tags_to_append: set[str] | None = None,
    replace: bool = False,
    existing_tags: dict[str, list] | None = None,
) -> dict[str, list]:
    """
    Process SAM tags into a mapping from tag name to [type, value].

    Args
    - sam_tags_list: list of SAM tags represented as tuples of (TAG, TYPE, VALUE)
    - tags_to_append: set of tag names for which to append values to existing tags
    - replace: if True, for tags not in tags_to_append, replace existing tags with the same name; if False, raise an
        error if a duplicate tag is found
    - existing_tags: existing mapping from tag name to [type, value] to update; if None, start with an empty mapping

    Returns: mapping from tag name to [type, value]
    """
    if tags_to_append is None:
        tags_to_append = set()
    if existing_tags is None:
        existing_tags = dict()
    for tag, value_type, value in sam_tags_list:
        if tag in existing_tags and tag in tags_to_append:
            value_type_existing, value_existing = existing_tags[tag]
            if value_type_existing == value_type and isinstance(value_existing, str) and isinstance(value, str):
                existing_tags[tag][1] = f"{value_existing}-{value}"
            else:
                raise ValueError(
                    f"Error: Cannot append ({value_type}, {value}) to existing tag {tag} with type, value of "
                    f"({value_type_existing}, {value_existing})."
                )
        elif replace or tag not in existing_tags:
            existing_tags[tag] = [value_type, value]
        else:
            raise ValueError(f"Duplicate SAM tag '{tag}' found.")
    return existing_tags


def read_parse(file: typing.IO, input_fmt: typing.Literal['fasta', 'fastq', 'bam']) -> \
    collections.abc.Generator[tuple[pysam.AlignedSegment | None, str, str, str | None], None, None]:
    """
    Parse reads from FASTA/FASTQ/BAM file.

    Args
    - file: file-like object to read from
    - input_fmt: 'fasta', 'fastq', or 'bam'

    Yields: read, read name, sequence, quality string
    - read is None for FASTA/FASTQ input; pysam.AlignedSegment for BAM input
    - quality string is None for FASTA input
    """
    if input_fmt == 'fasta':
        for name, seq in helpers.fasta_parse(file):
            yield None, name, seq, None
    elif input_fmt == 'fastq':
        for name, seq, _, qual in helpers.fastq_parse(file):
            yield None, name, seq, qual
    else:
        for read in file.fetch(until_eof=True):
            yield read, read.query_name, read.query_sequence, read.query_qualities_str


def process_barcode(
    barcode: list[str],
    num_tags: tuple[int, int, int],
    read_type_prefixes: set[str],
    reverse_barcode_order: bool,
) -> tuple[list[str], str | int, int]:
    """
    Process the barcode string from the read name into the final barcode and read type.

    Barcode tags in the read name are assumed to be ordered 5'->3' from each read orientation. For example, if 2 tags
    are sequenced in R1 and 4 tags are sequenced in R2, with 1 overlapping tag, then the barcode would look like
        [tag1][tag2][tag5][tag4][tag3][tag2]
    and num_tags should be (2, 4, 1). Overlapping tags are used to check for consistency. The final barcode is taken
    from the NUM_TAGS_IN_R1 tags from the R1 read and the first (NUM_TAGS_IN_R2 - NUM_TAGS_OVERLAP) tags from the R2
    read. If NUM_TAGS_TAGS_IN_R1 + NUM_TAGS_IN_R2 - NUM_TAGS_OVERLAP is less than the total number of tags in the
    barcode, then the extra tags at the end of the barcode are ignored.

    Args
    - barcode: list of tags in the barcode
    - num_tags: tuple of (NUM_TAGS_IN_R1, NUM_TAGS_IN_R2, NUM_TAGS_OVERLAP)
    - read_type_prefixes: set of prefixes to identify read type tags
    - reverse_barcode_order: whether to reverse the order of R2 barcode tags when adding to the SAM/BAM barcode tag

    Returns
    - barcode_final: final barcode as a list of tags
    - read_type: read type string; 0 if read_type_prefixes is an empty set, or -1 if no read type tag found despite
        prefixes
    - consistent: 1 if R1/R2 tags are consistent, -1 if inconsistent, or 0 if NUM_TAGS_OVERLAP == 0
    """
    n_r1, n_r2, n_overlap = num_tags
    R1_tags = barcode[:n_r1]
    R2_tags = barcode[n_r1:(n_r1 + n_r2)]
    if reverse_barcode_order:
        barcode_final = R1_tags + (R2_tags[:(n_r2 - n_overlap)])[::-1]
    else:
        barcode_final = R1_tags + R2_tags[:(n_r2 - n_overlap)]

    # determine R1-R2 tag consistency
    if n_overlap > 0:
        consistent = 1
        for i in range(n_overlap):
            if R1_tags[n_r1 - n_overlap + i] != R2_tags[n_r2 - i - 1]:
                # R1 and R2 barcodes do not match
                consistent = -1
                break
    else:
        # no overlapping tags to check
        consistent = 0

    # determine read type
    if len(read_type_prefixes) == 0:
        # no read type prefixes specified
        read_type = 0
    else:
        read_type = []
        for prefix in read_type_prefixes:
            barcode_final_new = []
            for tag in barcode_final:
                if tag.startswith(prefix):
                    read_type.append(tag)
                else:
                    barcode_final_new.append(tag)
            barcode_final = barcode_final_new
        if len(read_type) == 0:
            # no read type tags found
            read_type = -1
        else:
            read_type = '-'.join(read_type)

    return barcode_final, read_type, consistent


def convert_reads(
    file_in: str | typing.IO,
    file_out: str | typing.IO,
    file_discard: str | typing.IO | None,
    num_tags: tuple[int, int, int],
    input_fmt: typing.Literal['fasta', 'fastq', 'bam'],

    # read processing options
    remove_barcode_from_names: bool = True,
    reverse_barcode_order: bool = False,
    add_sample_to_barcode: str | None = None,
    add_SAM_tags: collections.abc.Iterable[str] | None = None,
    tag_read_type: str = "RT",
    tag_barcode: str = "CB",
    read_type_prefixes: collections.abc.Collection[str] | None = None,
    extract_UMI: tuple[int, int, bool] | None = None,

    # read filtering options
    discard_inconsistent_R1_R2: bool = False,
    require_read_type: bool = False,
    discard_UMI_N: bool = False,
    discard_UMI_mismatch: bool = False,

    # output options
    drop_quals: bool = False,
    uncompressed: bool = False,
    header: pysam.AlignmentHeader | dict | None = None,

    # performance options
    threads: int = 1,
    verbose: bool = False,
) -> dict[str, typing.Union[int, float]]:
    """
    Convert FASTA/FASTQ reads to BAM format with barcode and UMI tags.

    Args
    - file_in: input FASTA/FASTQ file path or file-like object
    - file_out: output BAM file path or file-like object
    - file_discard: output BAM file path or file-like object for discarded reads; if None, do not output discarded
    - num_tags: tuple of (NUM_TAGS_IN_R1, NUM_TAGS_IN_R2, NUM_TAGS_OVERLAP) specifying number of barcode tags from each
        read orientation
    - is_fasta: whether the input file is in FASTA format (as opposed to FASTQ)
    - remove_barcode_from_names: remove the barcode from read names after copying to the barcode and read type tags
    - reverse_barcode_order: reverse the order of R2 barcode tags when adding to the SAM/BAM barcode tag
    - add_sample_to_barcode: sample name to append to the barcode
    - add_SAM_tags: set of additional SAM tags to add to all reads, in TAG:TYPE:VALUE format
    - tag_read_type: 2-letter SAM/BAM tag name for the read type
    - tag_barcode: 2-letter SAM/BAM tag name for the barcode
    - read_type_prefixes: set of prefixes to identify read type tags in the barcode
    - extract_UMI: tuple of (start, stop, keep) for UMI extraction; if None, do not extract UMI
    - discard_inconsistent_R1_R2: (only relevant if num_tags[0] > 0) discard reads with inconsistent tags identified
        from R1 and R2 orientations
    - require_read_type: (only relevant if read_type_prefixes is not None) discard reads without a read type tag in the
        barcode
    - discard_UMI_N: discard reads with 'N' in the UMI sequence
    - discard_UMI_mismatch: discard reads where UMI sequences extracted from different sources do not match
    - drop_quals: (only for FASTQ input) do not store FASTQ quality scores in the BAM file
    - uncompressed: write uncompressed BAM output
    - header: BAM header to use; if None, a default header is created
    - threads: number of BAM output compression threads to use
    - verbose: print progress to stderr
    """
    # process arguments
    if header is None:
        header = dict(HD=dict(VN='1.6', SO='unknown'))
    barcode_suffix = f'.{add_sample_to_barcode}' if add_sample_to_barcode is not None else ''
    assert all(x >= 0 for x in num_tags), "Error: num_tags values must be non-negative."
    assert len(num_tags) == 3, "Error: num_tags must be a tuple of length 3."
    assert num_tags[2] <= min(num_tags[0], num_tags[1]), (
        "Error: NUM_TAGS_OVERLAP must be less than or equal to min(NUM_TAGS_IN_R1, NUM_TAGS_IN_R2)."
    )
    if read_type_prefixes is None:
        read_type_prefixes = set()
    elif not isinstance(read_type_prefixes, set):
        len_old = len(read_type_prefixes)
        read_type_prefixes = set(read_type_prefixes)
        if len(read_type_prefixes) < len_old and verbose:
            print(f"Warning: Duplicate entries found in read_type_prefixes: {read_type_prefixes}.", file=sys.stderr)
    regex_barcode = re.compile((num_tags[0] + num_tags[1]) * r"\[([a-zA-Z0-9_-]+)\]")
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1

    umi_start, umi_end, umi_keep = extract_UMI if extract_UMI is not None else (None, None, None)

    extra_SAM_tags_dict = process_SAM_tags(list(map(process_SAM_tag_str, add_SAM_tags))) if add_SAM_tags else dict()
    if 'N' in extra_SAM_tags_dict.get('RX', (None, ''))[1]:
        print("Warning: 'N' found in user-specified RX tag via --add_SAM_tags.", file=sys.stderr)

    if discard_UMI_mismatch:
        discard_UMI_N = True

    # initialize counts variables
    n_total = -1
    n_written = 0
    n_umi_mismatch = 0
    n_umi_mismatch = 0
    n_no_read_type = 0 if read_type_prefixes is not None else None
    n_r1_r2_mismatch = 0 if num_tags[2] > 0 else None

    with contextlib.ExitStack() as stack:
        f_in = stack.enter_context(open_reads_file(
            file_in,
            mode='rt' if input_fmt in ('fasta', 'fastq') else 'r',
            fmt=input_fmt,
            check_sq=False
        ))
        f_out = stack.enter_context(pysam.AlignmentFile(
            file_out,
            mode=mode_out,
            header=f_in.header if input_fmt == 'bam' else header,
            threads=threads
        ))
        f_discard = stack.enter_context(pysam.AlignmentFile(
            file_discard,
            mode=mode_out,
            header=f_in.header if input_fmt == 'bam' else header,
            threads=threads
        )) if file_discard is not None else None
        for n_total, (read, rname, seq, quals) in tqdm(enumerate(read_parse(f_in, input_fmt)), disable=not verbose):
            discard_reasons = set()
            quals = None if drop_quals else quals

            # extract UMI if requested
            umi_seq = None
            umi_qual = None
            UMI_tags = None
            if extract_UMI is not None:
                umi_seq = seq[umi_start:umi_end]
                if not f_discard and discard_UMI_N and 'N' in umi_seq:
                    n_umi_mismatch -= 1
                    continue
                UMI_tags = [('RX', 'Z', umi_seq)]
                if quals is not None:
                    umi_qual = quals[umi_start:umi_end]
                    UMI_tags.append(('QX', 'Z', umi_qual))
                if not umi_keep:
                    seq = seq[:umi_start] + seq[umi_end:]
                    if quals is not None:
                        quals = quals[:umi_start] + quals[umi_end:]

            # parse read name into components
            match = REGEX_READ_NAME.match(rname)
            if match is None:
                raise ValueError(f"Error: Invalid read name '{rname}'. Expected '<name>::<barcode>[ <SAM_tags>]'.")

            # process read name to output query name
            if remove_barcode_from_names:
                qname = match.group('name')
            else:
                qname = f"{rname}::{match.group('barcode')}"

            # process barcode
            if sum(num_tags) == 0:
                barcode_final = None
                read_type = 0 if len(read_type_prefixes) == 0 else -1
                consistent = 0
            else:
                barcode = list(regex_barcode.match(match.group('barcode')).groups())
                barcode_final, read_type, consistent = process_barcode(
                    barcode,
                    num_tags,
                    read_type_prefixes,
                    reverse_barcode_order
                )

            # apply read filtering on read type tags
            if read_type == -1:
                # no read type tag found, despite prefixes being specified
                if require_read_type:
                    n_no_read_type -= 1
                    if not f_discard:
                        continue
                    discard_reasons.add('no_read_type')
                else:
                    n_no_read_type += 1

            # apply read filtering on R1/R2 consistency
            if consistent == -1:
                # R1 and R2 barcodes do not match
                if discard_inconsistent_R1_R2:
                    n_r1_r2_mismatch -= 1
                    if not f_discard:
                        continue
                    discard_reasons.add('R1_R2_mismatch')
                else:
                    n_r1_r2_mismatch += 1

            # process SAM tags, in ascending order of precedence:
            # 1. from read name
            # 2. from UMI extraction
            # 3. from additional user-specified tags
            SAM_tags_str = match.group('sam_tags')
            final_SAM_tags = process_SAM_tags(list(map(process_SAM_tag_str, SAM_tags_str.split()))) if SAM_tags_str else dict()
            if UMI_tags:
                process_SAM_tags(UMI_tags, tags_to_append={'RX', 'QX'}, replace=False, existing_tags=final_SAM_tags)
            final_SAM_tags.update(extra_SAM_tags_dict)

            # apply read filtering on UMI
            if 'RX' in final_SAM_tags:
                if 'N' in final_SAM_tags['RX'][1]:
                    if discard_UMI_N:
                        n_umi_mismatch -= 1
                        if not f_discard:
                            continue
                        discard_reasons.add('UMI_N')
                    else:
                        n_umi_mismatch += 1
                else:
                    all_umis = final_SAM_tags['RX'][1].split('-')
                    if len(all_umis) > 1 and len(set(all_umis)) != len(all_umis):
                        # multiple UMIs extracted and do not match
                        if discard_UMI_mismatch:
                            n_umi_mismatch -= 1
                            if not f_discard:
                                continue
                            discard_reasons.add('UMI_mismatch')
                        else:
                            n_umi_mismatch += 1

            # create/modify read object
            if read is None:
                # convert FASTX read to SAM read
                read = pysam.AlignedSegment(f_out.header)
                read.flag = 4 # unmapped
                if len(seq) > 0:
                    read.query_sequence = seq
                if quals is not None and len(quals) > 0:
                    read.query_qualities_str = quals
            else:
                # modify existing SAM read
                if len(seq) != read.query_length:
                    read.query_sequence = seq
                if quals is not None and len(quals) != len(read.query_qualities_str):
                    read.query_qualities_str = quals
            read.query_name = qname

            if barcode_final is not None:
                barcode_final_str = '.'.join(barcode_final) + barcode_suffix
                read.set_tag(tag_barcode, barcode_final_str, value_type="Z")
            if read_type != 0:
                read.set_tag(tag_read_type, read_type, value_type="Z")
            for tag, (value_type, value) in final_SAM_tags.items():
                read.set_tag(tag, value, value_type=value_type if value_type != 'B' else None)
                # for array tags, pysam infers the type automatically and currently (v0.23.3) does not support passing
                # the value_type explicitly. See https://github.com/pysam-developers/pysam/issues/1108.

            if len(discard_reasons) > 0 and f_discard is not None:
                read.set_tag('XD', ';'.join(sorted(discard_reasons)), value_type='Z')
                f_discard.write(read)
            else:
                f_out.write(read)
                n_written += 1
        n_total += 1

    n_discard = sum([x for x in (n_umi_mismatch, n_no_read_type, n_r1_r2_mismatch) if x is not None and x < 0])
    assert n_total + n_discard == n_written, "Error: Read counts do not add up."

    counts = {
        'total': n_total,
        'written': n_written,
        'r1_r2_mismatch': n_r1_r2_mismatch,
        'no_read_type': n_no_read_type,
        'umi_mismatch': n_umi_mismatch,
    }
    return counts


def main():
    args = parse_arguments()

    if args.extract_UMI:
        umi_parts = args.extract_UMI.split(':')
        if len(umi_parts) != 3:
            raise ValueError("Error: Invalid --extract_UMI format. Expected start:stop:#.")
        umi_start = int(umi_parts[0])
        umi_end = int(umi_parts[1])
        umi_keep = bool(int(umi_parts[2]))
        extract_UMI = (umi_start, umi_end, umi_keep)
    else:
        extract_UMI = None

    counts = convert_reads(
        args.input if args.input else sys.stdin,
        args.output if args.output else sys.stdout,
        args.discarded,
        args.num_tags,
        input_fmt=args.input_fmt,

        remove_barcode_from_names=args.remove_barcode_from_names,
        reverse_barcode_order=args.reverse_barcode_order,
        add_sample_to_barcode=args.add_sample_to_barcode,
        add_SAM_tags=args.add_SAM_tags,
        tag_read_type=args.tag_read_type,
        tag_barcode=args.tag_barcode,
        read_type_prefixes=args.read_type_prefixes,
        extract_UMI=extract_UMI,

        # read filtering options
        discard_inconsistent_R1_R2=args.discard_inconsistent_R1_R2,
        require_read_type=args.require_read_type,
        discard_UMI_N=args.discard_UMI_N,
        discard_UMI_mismatch=args.discard_UMI_mismatch,

        # output options
        drop_quals=args.drop_quals,
        uncompressed=args.uncompressed,
        header=None,

        # performance options
        threads=args.threads,
        verbose=args.verbose,
    )

    if args.output_counts:
        with open(args.output_counts, 'w') as f:
            json.dump(counts, f)


if __name__ == "__main__":
    main()
