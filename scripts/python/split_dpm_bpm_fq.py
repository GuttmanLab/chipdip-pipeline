import gzip
import os
import argparse
import re
import pandas as pd
from collections import defaultdict
from helpers import fastq_parse, file_open

"""
Program to split barcoded reads into two files based on DPM and BPM tags and remove reads with impossible or incomplete barcods
"""


def parse_args():

    parser = argparse.ArgumentParser(
        description="Split fastq based on DPM or BPM barcode"
    )
    parser.add_argument(
        "--r1", dest="read_1", type=str, required=True, help="Fastq read 1"
    )
    parser.add_argument(
        "--format",
        dest="format",
        type=str,
        required=True,
        help="File with allowed barcodes for each position",
    )
    opts = parser.parse_args()

    return opts


def main():

    opts = parse_args()

    read_1_path = opts.read_1

    formatdict = load_format(opts.format)

    # Correctly formated DPM reads
    dpm_out_path = (
        os.path.splitext(os.path.splitext(read_1_path)[0])[0] + "_dpm.fastq.gz"
    )
    # Correctly formated BPM reads
    bpm_out_path = (
        os.path.splitext(os.path.splitext(read_1_path)[0])[0] + "_bpm.fastq.gz"
    )
    # Reads with barcode in incorrect order
    other_out_path = (
        os.path.splitext(os.path.splitext(read_1_path)[0])[0] + "_other.fastq.gz"
    )
    # Reads with NOT_FOUND barcode
    short_out_path = (
        os.path.splitext(os.path.splitext(read_1_path)[0])[0] + "_short.fastq.gz"
    )

    dpm_count = 0
    bpm_count = 0
    other_count = 0
    incomplete = 0
    counter = 0

    pattern = re.compile("\[([a-zA-Z0-9_\-]+)\]")

    with file_open(read_1_path) as read_1, \
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
                indexed = [i for i, t in enumerate(barcodes[1:]) if i != formatdict[t]]
                if len(indexed) != 0:
                    other_count += 1
                    other_out.write(
                        qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n"
                    )
                elif "DPM" in qname:
                    dpm_count += 1
                    dpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                elif "BEAD" in qname:
                    bpm_count += 1
                    bpm_out.write(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                else:
                    print("Read with unexpected barcode")
                    print(qname + "\n" + seq + "\n" + thrd + "\n" + qual + "\n")
                    raise Exception

    print("Reads without full barcode:", incomplete)
    print("DPM reads out:", dpm_count)
    print("BPM reads out:", bpm_count)
    print("Reads with incorrect barcode format:", other_count)


def load_format(formatfile):
    """
    Load file containing information on which barcodes can appear at which read positions
    """
    df = pd.read_csv(formatfile, header=None, sep="\t")
    return df.set_index(1)[0].to_dict(into=defaultdict(lambda: -1))


if __name__ == "__main__":
    main()
