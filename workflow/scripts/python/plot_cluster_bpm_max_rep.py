import argparse
import os

import pandas as pd
import matplotlib.figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# matplotlib note: here we bypass pyplot entirely to avoid having to choose a backend depending on whether this module
# is run interactively or not. See the following:
# - https://github.com/matplotlib/matplotlib/issues/26362
# - https://matplotlib.org/devdocs/gallery/user_interfaces/web_application_server_sgskip.html

"""
Generate maximum representation ecdfs for bead type representation within clusters.
"""


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the distributions of the counts and proportions of the maximally represented antibody ID over "
            "clusters to check for antibody ID uniqueness within clusters."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="<sample>.stats_bpm_max_rep.tsv.gz",
        help=(
            "Histogram of the BPM read count and proportion of the maximally represented antibody ID in each cluster, "
            "for each antibody ID. 4-column tab-delimited table with header: antibody ID, metric, value, and number "
            "of clusters. metric is 'count' or 'proportion', and value is the count or proportion of BPM reads from "
            "the maximally represented antibody ID in the cluster, with value rounded to 2 decimal places for "
            "proportions. Sample names are extracted from the filenames as the part before the first period, or "
            "supplied by the --samples option."
        )
    )
    parser.add_argument(
        "--sample_names",
        nargs="+",
        metavar="<sample_name>",
        help="Sample names to use for the input files. If used, must provide the same number of names as input files.",
    )
    parser.add_argument(
        "--counts",
        metavar="BPM_max_representation_counts.(png|pdf)",
        help="Path to plot of distribution of the read count of the maximally represented antibody ID in each cluster.",
    )
    parser.add_argument(
        "--proportions",
        metavar="BPM_max_representation_proportion.(png|pdf)",
        help="Path to plot of distribution of the proportion of the maximally represented antibody ID in each cluster.",
    )
    parser.add_argument(
        "-c",
        "--colorsequence",
        default="petroff10",
        help=(
            "Color sequence to use for different samples or read count categories. Supports options from matplotlib's "
            "color sequences (list(matplotlib.color_sequences.keys())) or anything that seaborn.color_palette() "
            "accepts, such as any matplotlib colormap (matplotlib.pyplot.colormaps()) or any built-in seaborn palette "
            "(deep, muted, bright, pastel, dark, colorblind)."
        )
    )
    return parser.parse_args()


def plot_ecdfs_proportions(df: pd.DataFrame, palette) -> matplotlib.figure.Figure:
    """
    Plot the distribution of the proportion of the maximally represented antibody ID in each cluster.

    Args
    - df: must contain the columns 'sample', 'value', and 'number of clusters'
    - palette: color palette to use for the plot

    Returns
    - fig: Left axes = eCDF over all clusters. Right axes = filtered for clusters with > 0 BPMs.
    """
    fig = matplotlib.figure.Figure(figsize=(10, 4), constrained_layout=True)
    axs = fig.subplots(1, 2, gridspec_kw=dict(wspace=0.1))
    sns.ecdfplot(
        data=df,
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[0],
    )
    sns.move_legend(axs[0], loc="upper left")
    axs[0].set_xlabel("Proportion of BPM reads from maximally represented antibody ID", fontsize=9)
    axs[0].set_ylabel("Proportion of clusters", fontsize=9)
    axs[0].set_title("All clusters")
    sns.ecdfplot(
        data=df.loc[df['value'] > 0],
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[1],
    )
    sns.move_legend(axs[1], loc="upper left")
    axs[1].set_xlabel("Proportion of BPM reads from maximally represented antibody ID", fontsize=9)
    axs[1].set_ylabel("Proportion of clusters", fontsize=9)
    axs[1].set_title("Filtered for clusters with > 0 BPMs")
    return fig


def plot_ecdfs_counts(df: pd.DataFrame, palette) -> tuple[matplotlib.figure.Figure, matplotlib.figure.Figure]:
    """
    Plot the distribution of the read counts of the maximally represented antibody ID in each cluster.

    Args
    - df: must contain the columns 'sample', 'value', and 'number of clusters'
    - palette: color palette to use for the plot

    Returns
    - fig1: eCDFs including clusters with 0 BPMs. Left axes = all clusters. Right axes = zoomed in.
    - fig2: eCDFs excluding clusters with 0 BPMs. Left axes = all clusters. Right axes = zoomed in.
    """
    fig1 = matplotlib.figure.Figure(figsize=(10, 4), constrained_layout=True)
    axs = fig1.subplots(1, 2, gridspec_kw=dict(wspace=0.1))
    sns.ecdfplot(
        data=df,
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[0],
    )
    sns.move_legend(axs[0], loc="lower right")
    axs[0].set_xlabel("BPM read count of maximally represented antibody ID", fontsize=9)
    axs[0].set_ylabel("Proportion of clusters", fontsize=9)
    sns.ecdfplot(
        data=df,
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[1],
    )
    sns.move_legend(axs[1], loc="lower right")
    axs[1].set_xlim(0, 30)
    axs[1].set_title("Zoomed in")
    axs[1].set_xlabel("BPM read count of maximally represented antibody ID", fontsize=9)
    axs[1].set_ylabel("Proportion of clusters", fontsize=9)
    fig1.suptitle("All clusters")

    fig2 = matplotlib.figure.Figure(figsize=(10, 4), constrained_layout=True)
    axs = fig2.subplots(1, 2, gridspec_kw=dict(wspace=0.1))
    sns.ecdfplot(
        data=df.loc[df['value'] > 0],
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[0],
    )
    sns.move_legend(axs[0], loc="lower right")
    axs[0].set_xlabel("BPM read count of maximally represented antibody ID", fontsize=9)
    axs[0].set_ylabel("Proportion of clusters", fontsize=9)

    sns.ecdfplot(
        data=df.loc[df['value'] > 0],
        x="value",
        weights="number of clusters",
        hue="sample",
        palette=palette,
        ax=axs[1],
    )
    sns.move_legend(axs[1], loc="lower right")
    axs[1].set_xlim(0, 30)
    axs[1].set_title("Zoomed in")
    axs[1].set_xlabel("BPM read count of maximally represented antibody ID", fontsize=9)
    axs[1].set_ylabel("Proportion of clusters", fontsize=9)
    fig2.suptitle("Filtered for clusters with > 0 BPMs")
    return fig1, fig2


def main():
    args = parse_arguments()
    assert args.counts or args.proportions, \
        "ERROR: At least one of the following arguments must be provided: --counts, --proportions"
    if args.sample_names is not None:
        if len(args.inputs) != len(args.sample_names):
            raise ValueError("ERROR: The number of input files must match the number of sample names.")
        sample_names = args.sample_names
    else:
        sample_names = [os.path.basename(f).split(".")[0] for f in args.inputs]
    if args.colorsequence in list(matplotlib.color_sequences.keys()):
        palette = sns.color_palette(matplotlib.color_sequences[args.colorsequence])
    else:
        palette = sns.color_palette(args.colorsequence)

    # read in cluster statistics from all samples
    dtype = {
        "antibody_id": 'category',
        "metric": 'category',
        "value": float,
        "number of clusters": int,
    }
    df = []
    for sample, path in zip(sample_names, args.inputs):
        df.append(
            pd.read_csv(path, sep="\t", header=0, dtype=dtype)
            .assign(sample=sample)
            .astype(dict(sample='category'))
        )
    df = pd.concat(df, axis=0, ignore_index=True)

    df_agg = (
        df.groupby(["sample", "metric", "value"], observed=True)
        ["number of clusters"].sum()
        .reset_index()
    )
    df_agg = pd.concat(
        (
            df_agg,
            (
                df_agg
                .groupby(["metric", "value"], observed=True)
                ["number of clusters"].sum()
                .reset_index()
                .assign(sample="all")
            )
        ),
        axis=0,
        ignore_index=True
    ).sort_values(["sample", "metric", "value"])

    if args.counts:
        fig1, fig2 = plot_ecdfs_counts(df_agg.loc[df_agg["metric"] == "count"], palette)
        with PdfPages(args.counts) as pdf:
            pdf.savefig(fig1, bbox_inches="tight")
            pdf.savefig(fig2, bbox_inches="tight")

    if args.proportions:
        fig = plot_ecdfs_proportions(df_agg.loc[df_agg["metric"] == "proportion"], palette)
        fig.savefig(args.proportions, bbox_inches="tight")


if __name__ == "__main__":
    main()
