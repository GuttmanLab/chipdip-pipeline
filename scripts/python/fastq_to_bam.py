import argparse
import os
import re
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import fastq_parse, file_open
import pandas as pd
import pysam

# regular expression for bead oligo name
PATTERN = re.compile("\[BEAD_([a-zA-Z0-9_\-]+)\]")


def main():
    args = parse_arguments()
    header = construct_sam_header(args.config)
    convert_reads(args.input, args.output, header, args.UMI_length)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert an antibody oligo FASTQ file to a BAM file")
    parser.add_argument("input", metavar="in.fastq", help="Input FASTQ file with barcodes in the read name")
    parser.add_argument("output", metavar="out.bam", help="Output BAM file")
    parser.add_argument("config", metavar="config.txt", help="Barcode config file")
    parser.add_argument("UMI_length", metavar="#", type=int, help="Length of the UMI")
    return parser.parse_args()


def construct_sam_header(config):
    """
    Args
    - config: str
        Path to a barcode config file

    Returns
    - header: pysam.AlignmentHeader
    """
    df = pd.read_csv(
        config,
        comment="#",
        on_bad_lines="warn",
        sep="\t",
        names=["Tag", "Name", "Sequence", "Number"],
    ).dropna()
    reference_names = [x.replace("BEAD_", "") for x in df["Name"] if "BEAD_" in x]
    header = pysam.AlignmentHeader().from_references(
        reference_names=reference_names,
        reference_lengths=[44444444] * len(reference_names),
    )
    return header


def convert_reads(path_in, path_out, header, UMI_length):
    """
    Args
    - path_in: str
        Path to input FASTQ file
    - path_out: str or file object
        Path to output BAM file, or a file object to write to
    - header: dict or pysam.AlignmentHeader
        SAM/BAM format header
    - UMI_length: int
        Length of the UMI
    """
    counter = 0
    with pysam.AlignmentFile(path_out, "wb", header=header) as output_bam, file_open(path_in) as reads:
        for qname, seq, _, _ in fastq_parse(reads):
            counter += 1
            if counter % 100000 == 0:
                print(counter, file=sys.stderr)
            match = PATTERN.search(qname)
            target_name = list(match.groups())[0]
            aligned_segment = initialize_alignment(header, qname, target_name, seq, UMI_length)
            output_bam.write(aligned_segment)
    print("The total number of bead reads was: ", counter, file=sys.stderr)


def initialize_alignment(header, query_name, reference_name, query_sequence, UMI_length):
    """
    Create a `pysam.AlignedSegment` object.

    The UMI is extracted from the first `UMI_length` bases, encoded as an integer using sequence_to_int(),
    and used as the reference start (column 4 in the SAM/BAM format).

    Args
    - header: pysam.AlignmentHeader
    - query_name: str
    - reference_name: str
    - query_sequence: str
    - UMI_length: int

    Returns
    - aligned_segment: pysam.AlignedSegment
    """
    aligned_segment = pysam.AlignedSegment(header)
    aligned_segment.query_name = query_name
    aligned_segment.reference_name = reference_name
    aligned_segment.reference_start = sequence_to_int(query_sequence[0:UMI_length], query_name)
    aligned_segment.query_sequence = query_sequence
    aligned_segment.flag = 0
    aligned_segment.cigar = ((0, len(query_sequence)),)
    return aligned_segment


def sequence_to_int(seq, query_name):
    table = "".maketrans({"A": "1", "C": "2", "G": "3", "T": "4", "N": "0"})
    try:
        value = int(seq.translate(table))
    except ValueError:
        print("The sequence ", query_name, "cannot be translated", file=sys.stderr)
        value = 0
    return value


if __name__ == "__main__":
    main()
