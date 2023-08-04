import argparse
import pysam
import pandas as pd
import re
from helpers import fastq_parse, file_open


def main():
    args = parse_arguments()
    header = construct_sam_header(args.config)
    convert_reads(args, header)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert a fastq file to a bam file")
    parser.add_argument(
        "-i", "--input", action="store", metavar="FILE", help="Input fastq"
    )
    parser.add_argument(
        "-o", "--output", action="store", metavar="FILE", help="Output BAM file"
    )
    parser.add_argument(
        "--config", metavar="CONFIG", action="store", help="The config file"
    )
    return parser.parse_args()


def construct_sam_header(config):
    df = pd.read_csv(
        config,
        comment="#",
        on_bad_lines="warn",
        sep="\t",
        names=["Tag", "Name", "Sequence", "Number"],
    ).dropna()
    names = list(df["Name"])
    beads = [x for x in names if "BEAD" in x]
    reference_names = list()
    for protein in beads:
        reference_names.append(protein.replace("BEAD_", ""))
    header = pysam.AlignmentHeader().from_references(
        reference_names=reference_names,
        reference_lengths=[44444444] * len(reference_names),
    )
    return header


def convert_reads(args, header):
    output_bam = pysam.AlignmentFile(args.output, "wb", header=header)
    pattern = re.compile("\[BEAD_([a-zA-Z0-9_\-]+)\]")
    counter = 0
    with file_open(args.input) as reads:
        for qname, seq, thrd, qual in fastq_parse(reads):
            counter = counter + 1
            if counter % 100000 == 0:
                print(counter)
            match = pattern.search(qname)
            protein_name = list(match.groups())[0]
            a = initialize_alignment(header, qname, protein_name, seq)
            output_bam.write(a)
    output_bam.close()
    print("The total number of bead reads was: ", counter)


def initialize_alignment(header, query_name, reference_name, query_sequence):
    """
    Create a `Pysam.AlignedSegment` object.

    The UMI is extracted from the first 8 bases, encoded as an integer using sequence_to_int(),
    and used as the reference start (column 4 in the SAM/BAM format).
    """
    a = pysam.AlignedSegment(header)
    a.query_name = query_name
    a.reference_name = reference_name
    a.reference_start = sequence_to_int(query_sequence[0:8], query_name)
    a.query_sequence = query_sequence
    a.flag = 0
    a.cigar = ((0, len(query_sequence)),)
    return a


def sequence_to_int(seq, query_name):
    table = "".maketrans({"A": "1", "C": "2", "G": "3", "T": "4", "N": "0"})
    try:
        value = int(seq.translate(table))
    except ValueError:
        print("The sequence ", query_name, "cannot be translated")
        value = 0
    return value


if __name__ == "__main__":
    main()
