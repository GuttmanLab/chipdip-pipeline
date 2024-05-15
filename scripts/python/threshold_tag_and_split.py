import argparse
import pysam
from collections import defaultdict, Counter
from pathlib import Path
import tqdm

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
        "-c",
        "--cluster_bam",
        dest="cluster_bam",
        type=str,
        required=True,
        help="Master aligned cluster Bamfile",
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
        "-d",
        "--dir",
        dest="dir",
        type=str,
        action="store",
        required=True,
        help="Directory to write individual antibody bams to",
    )
    parser.add_argument(
        "--min_oligos",
        action="store",
        type=int,
        required=True,
        help="The minimum number of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--proportion",
        action="store",
        type=float,
        required=True,
        help="The maximum representation proportion of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--max_size",
        action="store",
        type=int,
        required=True,
        help="The maximum cluster size to keep",
    )

    return parser.parse_args()

def main():
    args = parse_args()
    print("Using min_oligos: ", args.min_oligos)
    print("Using proportion: ", args.proportion)
    print("Using max_size: ", args.max_size)
    print("Writing bam to: ", args.output_bam)
    print("Writing splitbams to: ", args.dir)
    RG_dict = assign_labels(
        args.cluster_bam, args.min_oligos, args.proportion, args.max_size
    )
    label_bam_file(args.input_bam, args.output_bam, RG_dict)
    split_bam_by_RG(args.output_bam, args.dir)


def assign_labels(bamfile, min_oligos, threshold, max_size):

    # threshold_tag_and_split
    RG_dict = defaultdict(str)
    beads_dict = defaultdict(list)
    total_dict = defaultdict(int)

    with pysam.AlignmentFile(bamfile, 'rb') as inbam:
        
        for read in inbam.fetch(until_eof=True):
            read_type = read.get_tag("RT")
            barcode = read.get_tag("BC")
            chromesome = read.reference_name
            total_dict[barcode] += 1
            if "BEAD" in read_type or "BPM" in read_type:
                beads_dict[barcode].append(chromesome)

        # threshold tag and split
        for bc in total_dict.keys():
            count_beads = len(beads_dict[bc])
            if count_beads == 0:
                RG_dict[bc] = "none"
                continue
            cluster_size = total_dict[bc] - count_beads
            if int(cluster_size) > int(max_size):
                RG_dict[bc] = "filtered"
                continue
            bead_labels = Counter(beads_dict[bc])
            candidate = bead_labels.most_common()[0]
            if candidate[1] < min_oligos:
                RG_dict[bc] = "uncertain"
                continue
            elif candidate[1] / sum(bead_labels.values()) < threshold:
                RG_dict[bc] = "ambiguous"
                continue
            else:
                RG_dict[bc] = candidate[0]
                continue
            RG_dict[bc] = "malformed"

    return RG_dict

def label_bam_file(input_bam, output_bam, RG_dict):

    count, written, skipped = 0, 0, 0
    header = construct_read_group_header(input_bam, RG_dict)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for read in in_bam.fetch(until_eof=True):
            count += 1
            if count % 10000000 == 0:
                print(count)
            barcode = read.get_tag("BC")
            read_group = RG_dict[barcode]
            try:
                read.set_tag("RG", read_group, replace=True)
                out_bam.write(read)
                written += 1
            except KeyError:
                skipped += 1
    print("Total reads:", count)
    print("Reads written:", written)
    print("Reads with an error not written out:", skipped)

def construct_read_group_header(input_bam, RG_dict):
    """
    Add read group tags for each antibody to the bam file header

    Args:
        input_bam(str): Path to input master bam file
        RG_dict(dict): Dictionary of cluster assignments [barcode-> antibody]
    """
    proteins = set(RG_dict.values())
    sample_name = input_bam.split(".", 1)[0]
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        bam_header = input_file.header.to_dict()
        read_group_dict = [{"ID": name, "SM": sample_name} for name in list(proteins)]
        bam_header["RG"] = read_group_dict
        return bam_header
    
def split_bam_by_RG(output_bam, output_dir):
    """
    Split a bam file based on the read group

    Args:
        output_bam(str): Path to bam file containing reads with read group assignments
        output_dir(str): Directory to write read group specific bam files to
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    format_string = output_dir + "/%*_%!.%."
    pysam.split("-f", format_string, output_bam)

if __name__ == "__main__":
    main()