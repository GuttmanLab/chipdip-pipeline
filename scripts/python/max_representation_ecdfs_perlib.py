import argparse
from collections import Counter
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")

"""
Generate maximum representation ecdfs for bead type representation within clusters.
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    ecdf_plot_ax = "None"
    for f in files:
        ecdf_plot_ax = max_representation_ecdf(f, ecdf_plot_ax)
    ecdf_plot_fig = ecdf_plot_ax.get_figure()
    ecdf_plot_fig.savefig(
        args.directory + "/Max_representation_ecdf.pdf", bbox_inches="tight"
    )
    ecdf_counts_ax = "None"
    for f in files:
        ecdf_counts_ax = max_representation_ecdf_counts(f, ecdf_counts_ax, args.xlim)
    ecdf_counts_fig = ecdf_counts_ax.get_figure()
    ecdf_counts_fig.savefig(
        args.directory + "/Max_representation_counts.pdf", bbox_inches="tight"
    )


def max_representation_ecdf(clusterfile, ax):
    """
    Plot maximum representation ecdf for bead representation within clusters of a single clusterfile

    Args:
        clusterfile(str): path to clusterfile
        ax(obj): axis to plot on
    """
    results = []
    clustername = clusterfile.replace(".clusters", "")
    with open(clusterfile, "r") as clusters:
        for line in clusters:
            barcode, *reads = line.rstrip("\n").split("\t")
            bead_reads = [read for read in reads if read.startswith("BPM")]
            if len(bead_reads) > 1:  # ignore singleton clusters
                bead_labels = Counter(
                    [read.split(":")[0].split("_", 1)[1] for read in bead_reads]
                )
                candidate = bead_labels.most_common()[0]
                results.append(candidate[1] / len(bead_reads))
    if ax == "None":
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth=3, ax=ax, label=clustername)
    ax.set(
        xlabel="Maximum Bead Representation Proportion", ylabel="Proportion of Beads"
    )
    ax.legend()
    return ax


def max_representation_ecdf_counts(clusterfile, ax, xlimit):
    """
    Plot counts of BPM reads with maximum representation

    Args:
        clusterfile(str): path to clusterfile
        ax(obj): axis to plot on
        xlimit(int): maximum x value to show in plot
    """
    results = []
    clustername = clusterfile.replace(".clusters", "")
    with open(clusterfile, "r") as clusters:
        for line in clusters:
            barcode, *reads = line.rstrip("\n").split("\t")
            bead_reads = [read for read in reads if read.startswith("BPM")]
            if len(bead_reads) > 1:  # Exclude singleton clusters
                bead_labels = Counter(
                    [read.split(":")[0].split("_", 1)[1] for read in bead_reads]
                )
                candidate = bead_labels.most_common()[0]
                results.append(candidate[1])
    if ax == "None":
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth=3, ax=ax, label=clustername)
    ax.set(xlabel="Num Oligos for Most Common Type", ylabel="Proportion of Beads")
    ax.set(xlim=(0, int(xlimit)))
    ax.legend()
    return ax


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the maximum representation ecdf plots to check for bead type uniqueness within clusters."
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        required=True,
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern",
        action="store",
        required=True,
        help="The pattern of cluster file names",
    )
    parser.add_argument(
        "--xlim",
        action="store",
        required=True,
        help="The maximum x value on the counts plot",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
