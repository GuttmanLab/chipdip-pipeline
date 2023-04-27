import tqdm
import glob
import argparse

"""
Count 1) number of clusters, 2) number of DPM reads (aligned) and 3) number of  BPM reads within each clusterfile for a directory of clusterfiles.
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    for f in files:
        count_statistics(f)


def count_statistics(clusterfile):
    """
    Loop through all clusters within a clusterfile, counting DPM and BPM reads

    Args:
        clusterfile(str): Path to clusterfile
    """
    cluster = 0
    dpm = 0
    bpm = 0
    with open(clusterfile, "r") as clusters:
        for line in tqdm.tqdm(clusters):
            cluster += 1
            dpm += line.count("DPM[")
            bpm += line.count("BPM[")
    print("For clusterfile ", clusterfile)
    print("Total number of clusters: ", cluster)
    print("Total number of BPM: ", bpm)
    print("Total number of DPM: ", dpm)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the statistics for all clusterfiles in a directory"
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern", action="store", help="The pattern of cluster file names"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
