import argparse
import pysam
import re
from collections import defaultdict
import os

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
    parser.add_argument(
        "-o",
        "--output_bam",
        dest="output_bam",
        type=str,
        required=True,
        help="Path to output master bam with antibody tag added",
    )
    parser.add_argument(
        "--num_tags",
        dest="num_tags",
        type=int,
        required=True,
        help="Number of tags in barcode",
    )

    return parser.parse_args()

def main():
    args = parse_args()
    print("Writing tagged bam to: ", args.output_bam)
    label_bam_file(args.input_bam, args.output_bam, args.num_tags)

def label_bam_file(input_bam, output_bam, num_tags):
    """
    Add antibody label to individual reads of the master DNA bam file

    Args:
        input_bam(str): Path to input master bam file
        output_bam(str): Path to write labeled bam file
        num_tags(int): number of tags in barcode
    """
    count, written, duplicates, skipped = 0, 0, 0, 0
    pattern = re.compile("::" + num_tags * "\[([a-zA-Z0-9_\-]+)\]")
    rt_pattern = re.compile(r"RPM|BPM|DPM|BEAD")
    found = defaultdict(set)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:
        for read in in_bam.fetch(until_eof=True):
            count += 1
            if count % 100000 == 0:
                print(count)
            name = read.query_name
            match = pattern.search(name)
            full_barcode = list(match.groups())
            read_type = rt_pattern.findall(full_barcode[0])[0]
            barcode = full_barcode[1:]
            ref_barcode = ".".join(barcode)
            if "DPM" in read_type:
                position = read.reference_name + ":" + str(read.reference_start) + '-' + str(read.reference_end)
            elif "BPM" in read_type or "BEAD" in read_type:
                position = read.reference_name + ":" + str(read.reference_start) + '-' + str(0)
            if position in found[ref_barcode]:
                duplicates += 1
            else:
                try:
                    found[ref_barcode].add(position)
                    read.set_tag("RT", read_type, replace=True)
                    read.set_tag("RC", ref_barcode, replace=True)
                    out_bam.write(read)
                    written += 1
                except KeyError:
                    skipped += 1

    # sort bam file by barcode
    # pysam.sort("-t", "BC", "-o", output_bam, output_bam)
    print("Total reads:", count)
    print("Reads written:", written)
    print("Duplicate reads:", duplicates)
    print("Reads with an error not written out:", skipped)

if __name__ == "__main__":
    main()