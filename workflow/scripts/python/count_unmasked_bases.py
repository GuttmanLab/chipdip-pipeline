"""
Count the number of non-masked nucleotides in a FASTA file.
"""
import argparse
from collections import Counter
import collections.abc
import os
import sys
import typing
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import file_open


def count_unmasked(file: typing.TextIO, unmasked_alphabet: collections.abc.Iterable[str] = ('A', 'C', 'G', 'T')) -> int:
    """
    Count the number of non-masked nucleotides in a FASTA file.

    All characters (except terminal whitespace) in lines not starting with '>' are considered.

    Args:
    - file: FASTA file
    - unmasked_alphabet: Letters to include in the unmasked alphabet. Case-sensitive.
    """
    c = Counter()
    for line in file:
        if not line.startswith('>'):
            c.update(Counter(line.rstrip()))
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
        help=(
            "Letters (case-sensitive) to include in the unmasked alphabet. "
            "Use many times (e.g., -u A -u C -u T -u G) to specify multiple letters. "
            "If not specified, defaults to A, C, G, and T."
        )
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    unmasked_alphabet = ('A', 'C', 'G', 'T') if args.unmasked_alphabet is None else args.unmasked_alphabet
    if args.input == '-':
        effective_genome_size = count_unmasked(sys.stdin, unmasked_alphabet)
    else:
        with file_open(args.input, mode="rt") as f:
            effective_genome_size = count_unmasked(f, unmasked_alphabet)
    print(effective_genome_size)
    return effective_genome_size


if __name__ == "__main__":
    main()
