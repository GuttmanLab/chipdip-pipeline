import argparse
import gzip
import os
import re
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import fastq_parse, file_open
import pandas as pd

"""
Program to split barcoded reads into two files based on DPM and BPM tags and remove reads with impossible or incomplete barcods
"""


def parse_args():

    parser = argparse.ArgumentParser(
        description="Split fastq based on DPM or BPM barcode"
    )
    parser.add_argument(
        "--r1", dest="read_1", type=str, required=True, help="Path to FASTQ file (optionally gzipped) for read 1"
    )
    parser.add_argument(
        "--format",
        metavar="FILE",
        type=str,
        help="(Recommended) Path to format file specifying allowed tags for each position in the barcode."
    )
    opts = parser.parse_args()

    return opts


def main():

    opts = parse_args()

    read_1_path = opts.read_1

    if opts.format:
        valid_tag_positions = load_format(opts.format)

    base_path = os.path.splitext(os.path.splitext(read_1_path)[0])[0]
    # Correctly formated DPM reads
    dpm_out_path = base_path + "_dpm.fastq.gz"
    # Correctly formated BPM reads
    bpm_out_path = base_path + "_bpm.fastq.gz"
    # Reads with barcode in incorrect order
    other_out_path = base_path + "_other.fastq.gz"
    # Reads with NOT_FOUND barcode
    short_out_path = base_path + "_short.fastq.gz"

    dpm_count = 0
    bpm_count = 0
    other_count = 0
    incomplete = 0
    counter = 0

    pattern = re.compile(r"\[([a-zA-Z0-9_\-]+)\]")

    with file_open(read_1_path, mode="rt") as read_1, \
         gzip.open(dpm_out_path, "wt") as dpm_out, \
         gzip.open(bpm_out_path, "wt") as bpm_out, \
         gzip.open(other_out_path, "wt") as other_out, \
         gzip.open(short_out_path, "wt") as short_out:

        for qname, seq, thrd, qual in fastq_parse(read_1):
            counter += 1
            barcodes = pattern.findall(qname.split('::')[1])
            if counter % 10000 == 0:
                print(counter)
            if "NOT_FOUND" in barcodes:
                incomplete += 1
                short_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
            else:
                unexpected_tag = False
                if opts.format:
                    for i, tag in enumerate(barcodes[1:]):
                        # first tag in barcodes should either be a DPM or BEAD tag;
                        # verify that subsequent tags are in valid positions
                        allowed_indices = valid_tag_positions.get(tag)
                        if allowed_indices is None:
                            print(f'Tag {tag} not in format file.', file=sys.stderr)
                            unexpected_tag = True
                            break
                        elif i not in allowed_indices:
                            unexpected_tag = True
                            break
                if unexpected_tag:
                    other_count += 1
                    other_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                elif "DPM" in qname:
                    dpm_count += 1
                    dpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                elif "BEAD" in qname:
                    bpm_count += 1
                    bpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                else:
                    print("Read with unexpected barcode", file=sys.stderr)
                    print(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n",
                          file=sys.stderr)
                    raise Exception

    print("Reads without full barcode:", incomplete)
    print("DPM reads out:", dpm_count)
    print("BPM reads out:", bpm_count)
    print("Reads with incorrect barcode format:", other_count)


def load_format(formatfile: str) -> dict[str, tuple[int]]:
    """
    Load file containing information on which tags can appear at which read positions

    Returns: mapping from tag name to expected read positions (rounds)
    """
    df = pd.read_csv(
        formatfile,
        sep="\t",
        header=None,
        names=["round", "name"],
        usecols=[0, 1],
        index_col=False,
    )
    d = (
        df.groupby("name")["round"]
        .apply(lambda s: tuple(sorted(s.unique())))
        .to_dict()
    )
    for name, rounds in d.items():
        if len(rounds) > 1 and -1 in rounds:
            print(
                (
                    f"The format file indicates tag {name} as both being used and "
                    f"not used (round = -1) in the experiment: {rounds}."
                ),
                file=sys.stderr,
            )
    return d


if __name__ == "__main__":
    main()
