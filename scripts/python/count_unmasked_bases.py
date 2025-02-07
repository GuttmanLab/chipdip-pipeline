"""
Count the number of non-masked nucleotides in a FASTA file.
"""

from collections import Counter
import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import file_open


def count_unmasked(file, unmasked_alphabet=None) -> int:
    """
    Count the number of non-masked nucleotides in a FASTA file.

    All characters (except terminal whitespace) in lines not starting with '>' are considered.

    Args:
    - file: file object (io.TextIOBase)
        FASTA file. Must be opened for reading in text mode.
    - unmasked_alphabet: Container. default=None
        Letters to include in the unmasked alphabet. Case-sensitive.
        If None, default to ['A', 'C', 'G', 'T'].
    """
    c = Counter()
    for line in file:
        if not line.startswith('>'):
            c.update(Counter(line.rstrip()))
    if unmasked_alphabet is None:
        unmasked_alphabet = ['A', 'C', 'G', 'T']
    return sum(c[nt] for nt in unmasked_alphabet)


def parse_arguments():
    parser = argparse.ArgumentParser(description=(
        "Count the number of non-masked nucleotides in a FASTA file. "
        "Writes calculated value to standard out."
    ))
    parser.add_argument(
        "input",
        metavar="in.fa(.gz)|-",
        help="Input FASTA file; use '-' for standard input."
    )
    parser.add_argument(
        "-u", "--unmasked_alphabet",
        metavar="A|C|T|G",
        action='append',
        default=[],
        help=(
            "Letters (case-sensitive) to include in the unmasked alphabet. "
            "Use many times (e.g., -u A -u C -u T -u G) to specify multiple letters. "
            "If not specified, defaults to A, C, G, and T."
        )
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if len(args.unmasked_alphabet) == 0:
        args.unmasked_alphabet = None
    if args.input == '-':
        effective_genome_size = count_unmasked(sys.stdin, args.unmasked_alphabet)
    else:
        with file_open(args.input, mode="rt") as f:
            effective_genome_size = count_unmasked(f, args.unmasked_alphabet)
    print(effective_genome_size)
    return effective_genome_size


if __name__ == "__main__":
    main()
