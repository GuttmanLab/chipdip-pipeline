import tqdm
import glob
import argparse
import pysam
import numpy as np
from collections import Counter, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import os
import pandas as pd

"""
Gnerate all cluster statistics from master bam file
Merging the function of generate_cluster_statistics.py, max_representation_ecdfs_perlib.py, get_bead_size_distribution.py

Specifically,
    Count 1) number of clusters, 2) number of DPM reads (aligned) and 3) number of  BPM reads within each clusterfile for a directory of clusterfiles.
    Generate maximum representation ecdfs for bead type representation within clusters.
    Profiles the proportion of clusters and proportion of reads within various cluster size categories. Considers DPM and BPM reads seperately (based on the input parameter readtype).
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    
    # 
    cluster_counts_dpm = []
    read_counts_dpm = []

    cluster_counts_bpm = []
    read_counts_bpm = []

    # max_representation_ecdf
    ecdf_plot_ax = "None"
    ecdf_counts_ax = "None"

    for f in files:
        ecdf_plot_ax, ecdf_counts_ax, df1_dpm, df2_dpm, df1_bpm, df2_bpm = generate_statistics(f, 
            ecdf_plot_ax, ecdf_counts_ax, args.xlim, args.directory)
        cluster_counts_dpm.append(df1_dpm)
        read_counts_dpm.append(df2_dpm)
        cluster_counts_bpm.append(df1_bpm)
        read_counts_bpm.append(df2_bpm)

    ecdf_plot_fig = ecdf_plot_ax.get_figure()
    ecdf_plot_fig.savefig(
        args.directory + "/Max_representation_ecdf.pdf", bbox_inches="tight"
    )
    ecdf_counts_fig = ecdf_counts_ax.get_figure()
    ecdf_counts_fig.savefig(
        args.directory + "/Max_representation_counts.pdf", bbox_inches="tight"
    )

    cluster_df_dpm = pd.concat(cluster_counts_dpm, axis=1).transpose()
    read_df_dpm = pd.concat(read_counts_dpm, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df_dpm)
    cluster_fig.savefig(
        args.directory + "/" + "DPM" + "_cluster_distribution.pdf",
        bbox_inches="tight",
    )
    read_fig = plot_profile(read_df_dpm)
    read_fig.savefig(
        args.directory + "/" + "DPM" + "_read_distribution.pdf",
        bbox_inches="tight",
    )

    cluster_df_bpm = pd.concat(cluster_counts_bpm, axis=1).transpose()
    read_df_bpm = pd.concat(read_counts_bpm, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df_bpm)
    cluster_fig.savefig(
        args.directory + "/" + "BPM" + "_cluster_distribution.pdf",
        bbox_inches="tight",
    )
    read_fig = plot_profile(read_df_bpm)
    read_fig.savefig(
        args.directory + "/" + "BPM" + "_read_distribution.pdf",
        bbox_inches="tight",
    )

def generate_statistics(bamfile, ax1, ax2, xlimit, dir):
    """
    Loop through all clusters within a bamfile, counting DPM and BPM reads

    Args:
        bamfile(str): Path to bamfile
    """
    # count statistics
    cluster = 0
    dpm = 0
    bpm = 0
    barcodes = set()

    # max_representation_ecdf
    results = []
    results_counts = []
    bamname = os.path.basename(bamfile)

    # get bead size distribution
    # count bpm and dpm at the same time
    bins = np.array([0, 1, 5, 10, 20, 50, 100, 200])  # bins with data are 1-8
    cluster_counts_dpm = defaultdict(int)
    read_counts_dpm = defaultdict(int)
    cluster_counts_bpm = defaultdict(int)
    read_counts_bpm = defaultdict(int)

    # the cluster dictionaries
    beads_dict = defaultdict(set)
    dpm_dict = defaultdict(int)

    count = 0
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            count += 1
            if count % 100000 == 0:
                print(count)
            read_type = read.get_tag("RT")
            barcode = read.get_tag("RC")
            try:
                chromesome = read.reference_name
                if "DPM" in read_type:
                    dpm_dict[barcode] += 1 # get bead size distribution
                elif "BEAD" in read_type or "BPM" in read_type:
                    beads_dict[barcode].add((chromesome, str(count)))
            except KeyError:
                pass

    for bin in np.arange(
        1, len(bins) + 1
    ):  # initialize all bins in case any end up being empty categories
        cluster_counts_dpm[bin] = 0
        read_counts_dpm[bin] = 0
        cluster_counts_bpm[bin] = 0
        read_counts_bpm[bin] = 0
    
    for bc in beads_dict.keys():
        beads_list = [ele[0] for ele in beads_dict[bc]]
        count_beads = len(beads_list)
        # max representation ecdf
        if count_beads > 1:
            candidate = Counter(beads_list).most_common()[0]
            results.append(candidate[1] / count_beads)
            results_counts.append(candidate[1])
        # get bead size distribution
        cluster += 1
        barcodes.add(bc)
        bpm += count_beads
        bin_bpm = np.digitize(count_beads, bins, right=True)
        cluster_counts_bpm[bin_bpm] += 1
        read_counts_bpm[bin_bpm] += count_beads
    
    # get bead size distribution
    for bc in dpm_dict.keys():
        count_dpm = dpm_dict[bc]
        if bc not in barcodes:
            cluster += 1
        dpm += count_dpm
        bin_dpm = np.digitize(count_dpm, bins, right=True)
        cluster_counts_dpm[bin_dpm] += 1
        read_counts_dpm[bin_dpm] += count_dpm


    # generate cluster statistics
    print("For bamfile ", bamfile)
    print("Total number of clusters: ", cluster)
    print("Total number of BPM: ", bpm)
    print("Total number of DPM: ", dpm)

    filepath = os.path.join(dir, "cluster_statistics.txt")
    f = open(filepath, "a")
    f.write("For bamfile " + str(bamfile) + "\n")
    f.write("Total number of clusters: " + str(cluster) + "\n")
    f.write("Total number of BPM: " + str(bpm) + "\n")
    f.write("Total number of DPM: " + str(dpm) + "\n")
    f.close()

    # max representation ecdf
    if ax1 == "None":
        ax1 = plt.figure().subplots()
    ax1 = sns.ecdfplot(results, linewidth=3, ax=ax1, label=bamname)
    ax1.set(
        xlabel="Maximum Bead Representation Proportion", ylabel="Proportion of Beads"
    )
    ax1.legend()
    if ax2 == "None":
        ax2 = plt.figure().subplots()
    ax2 = sns.ecdfplot(results_counts, linewidth=3, ax=ax2, label=bamname)
    ax2.set(xlabel="Num Oligos for Most Common Type", ylabel="Proportion of Beads")
    ax2.set(xlim=(0, int(xlimit)))
    ax2.legend()

    # get bead size distribution
    df_cluster_counts_dpm = pd.DataFrame.from_dict(cluster_counts_dpm, orient="index")
    df_cluster_counts_dpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]
    df_read_counts_dpm = pd.DataFrame.from_dict(read_counts_dpm, orient="index")
    df_read_counts_dpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]

    df_cluster_counts_bpm = pd.DataFrame.from_dict(cluster_counts_bpm, orient="index")
    df_cluster_counts_bpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]
    df_read_counts_bpm = pd.DataFrame.from_dict(read_counts_bpm, orient="index")
    df_read_counts_bpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]

    return ax1, ax2, df_cluster_counts_dpm, df_read_counts_dpm, df_cluster_counts_bpm, df_read_counts_bpm

def plot_profile(df):
    """
    Plot the proportion of reads within each size category as a stacked bar graph
    Args:
        df(dataframe): binned counts
    """
    columns = ["1", "2-5", "6-10", "11-20", "21-50", "51-100", "101-200", "201+"]
    df.columns = columns
    df = df.div(df.sum(axis=1), axis=0)
    plot = df.plot(
        kind="bar", stacked=True, ylabel="Proportion", cmap=cm.get_cmap("Dark2")
    )
    plot.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    plot.set_ylabel = "Proportion"
    return plot.get_figure()

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
    parser.add_argument(
        "--xlim",
        action="store",
        required=True,
        help="The maximum x value on the counts plot",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()