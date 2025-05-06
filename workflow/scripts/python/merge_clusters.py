import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import cluster as c

"""
Merges clusters containing the same barcode. Input cluster file is assumed to be in sorted order.
"""


def main():
    args = parse_arguments()
    c.merge_clusters(args.input, args.output)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merges clusters in a clusterfile that contain the same barcode"
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="FILE",
        action="store",
        required=True,
        help="The input clusters file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        action="store",
        required=True,
        help="The output clusters file.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
