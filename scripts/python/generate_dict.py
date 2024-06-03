import argparse
import pysam
import re
from collections import defaultdict
import json
import os
import pickle

'''
Python code to add read type and barcode tags to input bam file
'''

def parse_args():

    parser = argparse.ArgumentParser(
        description="Add antibody label to DNA bamfile and generate individual bamfiles for each antibody"
    )
    parser.add_argument(
        "-i",
        "--input_bam",
        dest="input_bam",
        type=str,
        required=True,
        help="Master aligned DNA Bamfile",
    )

    return parser.parse_args()

def main():
    args = parse_args()
    generate_dictionary(args.input_bam)

def generate_dictionary(input_bam):
    """
    Add antibody label to individual reads of the master DNA bam file

    Args:
        input_bam(str): Path to input master bam file
        output_bam(str): Path to write labeled bam file
        num_tags(int): number of tags in barcode
    """
    count = 0
    beads_dict = defaultdict(set)
    dpm_dict = defaultdict(int)
    total_dict = defaultdict(int)

    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        for read in in_bam.fetch(until_eof=True):
            count += 1
            if count % 1000000 == 0:
                print(count)
            read_type = read.get_tag("RT")
            barcode = read.get_tag("RC")
            try:
                chromesome = read.reference_name
                total_dict[barcode] += 1
                if "DPM" in read_type:
                    dpm_dict[barcode] += 1 # get bead size distribution
                elif "BEAD" in read_type or "BPM" in read_type:
                    beads_dict[barcode].add((chromesome, str(count)))
            except KeyError:
                pass
    
    final_dict = {"beads_dict": beads_dict,
                  "dpm_dict": dpm_dict,
                  "total_dict": total_dict}
    
    dictname = input_bam + ".dict"

    with open(dictname, "wb") as f:
        pickle.dump(final_dict, f)

if __name__ == "__main__":
    main()