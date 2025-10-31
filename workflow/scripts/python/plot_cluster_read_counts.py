import argparse
import collections.abc
import os

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.use("Agg")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot distributions of reads and clusters by reads per cluster for each read type."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="<sample>.stats_reads_per_cluster.tsv.gz",
        help=(
            "Histogram of reads per cluster, for each antibody and read type. 4-column tab-delimited table with "
            "header: antibody ID, read type, reads per cluster, number of clusters. Sample names are extracted from "
            "the filenames as the part before the first period, or supplied by the --samples option."
        )
    )
    parser.add_argument(
        "--sample_names",
        nargs="+",
        metavar="<sample_name>",
        help="Sample names to use for the input files. If used, must provide the same number of names as input files.",
    )
    parser.add_argument(
        "--dpm_read",
        metavar="DPM_read_distribution.(png|pdf)",
        help="Path to plot of distribution of DPM reads by the size (in terms of DPM reads) of their cluster.",
    )
    parser.add_argument(
        "--dpm_cluster",
        metavar="DPM_cluster_distribution.(png|pdf)",
        help="Path to plot of distribution of clusters by number of DPM reads.",
    )
    parser.add_argument(
        "--bpm_read",
        metavar="BPM_read_distribution.(png|pdf)",
        help="Path to plot of distribution of BPM reads by the size (in terms of BPM reads) of their cluster.",
    )
    parser.add_argument(
        "--bpm_cluster",
        metavar="BPM_cluster_distribution.(png|pdf)",
        help="Path to plot of distribution of clusters by number of BPM reads.",
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


# Plots
# - Total
#   - eCDF, no limits on x-axis, hue = sample
#   - eCDF, x-axis limited to 0-30, hue = sample
#   - binned barplots, x-axis = sample
#
# - Per antibody ID
#   - y-axis = antibody ID
#   - x-axis = reads per cluster, errorbars indicate 0.025 and 0.975 quantiles
#   - hue = sample


# antibody ID, read type, reads per cluster, number of clusters
# read type, reads per cluster, number of clusters


def cluster_distribution_ecdf(
    df: pd.DataFrame,
    read_type: str,
    variable: str,
    zoom: int | None = 30,
    palette: str = "petroff10",
) -> matplotlib.figure.Figure:
    """
    Plot the distribution of clusters or reads by number of reads per cluster.

    Args
    - df: DataFrame with columns "sample", "reads per cluster", and "number of clusters" or "number of reads".
        Should already be subsetted for the read type to plot.
    - read_type: The read type to plot (e.g., DPM, BPM).
    - variable: One of "clusters" or "reads" - whether to plot a distribution of clusters or reads.
    - zoom: The x-axis limit for the zoomed in plot. If None, no zoom is applied.
    - palette: The color palette to use for samples.

    Returns: figure object
    """
    if zoom:
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), gridspec_kw=dict(wspace=0.1), constrained_layout=True)
    else:
        fig, ax = plt.subplots(figsize=(5, 3.5), constrained_layout=True)
        axs = [ax]
    sns.ecdfplot(
        data=df,
        x="reads per cluster",
        weights=f"number of {variable}",
        hue='sample',
        palette=palette,
        ax=axs[0],
    )
    sns.move_legend(axs[0], "upper left", bbox_to_anchor=(1, 1))
    axs[0].set_xlabel(f"Number of {read_type} reads per cluster")
    label = variable if variable == "clusters" else f"{read_type} reads"
    title = f"Distribution (eCDF) of {label} by cluster size"
    axs[0].set_title(title, fontsize=10)

    if zoom:
        sns.ecdfplot(
            data=df,
            x="reads per cluster",
            weights=f"number of {variable}",
            hue='sample',
            palette=palette,
            ax=axs[1],
        )
        sns.move_legend(axs[1], "upper left", bbox_to_anchor=(1, 1))
        axs[1].set_xlim(-1, zoom)
        axs[1].set_xlabel(f"Number of {read_type} reads per cluster")
        axs[1].set_title(title + " (zoomed in)", fontsize=10)
    return fig


def cluster_distribution_stacked_bar(
    df: pd.DataFrame,
    read_type: str,
    variable: str,
    palette: collections.abc.Sequence[str] = matplotlib.color_sequences['petroff10'],
) -> matplotlib.figure.Figure:
    """
    Plot the distribution of clusters or reads by number of reads per cluster.

    Args
    - df: DataFrame with columns "sample", "reads per cluster", and "number of clusters" or "number of reads".
        Should already be subsetted for the read type to plot.
    - read_type: The read type to plot (e.g., DPM, BPM).
    - variable: One of "clusters" or "reads" - whether to plot a distribution of clusters or reads.
    - palette: The color palette to use for reads per cluster bins.

    Returns: figure object
    """
    col_weights = f"number of {variable}"
    label = variable if variable == "clusters" else f"{read_type} reads"
    df_wide = (
        df
        .groupby(["sample", "reads per cluster, binned"], observed=False)
        [col_weights].sum()
        .reset_index()
        # 3 relevant columns: sample, reads per cluster, number of clusters/reads
        .groupby("sample", observed=True)[["reads per cluster, binned", col_weights]]
        .apply(lambda g: g.assign(proportion=g[col_weights] / g[col_weights].sum()))
        .droplevel(1, axis=0)
        # 3 relevant columns: sample, reads per cluster, proportion
        .reset_index()
        .pivot(index="sample", columns="reads per cluster, binned", values="proportion")
    )
    n_samples = len(df_wide)
    fig, ax = plt.subplots(figsize=(n_samples + 2, 4), constrained_layout=True)
    x = np.arange(n_samples)
    bottoms = np.zeros(n_samples)
    for i, bin in enumerate(df_wide.columns):
        proportions = df_wide[bin].to_numpy()
        ax.bar(x, proportions, bottom=bottoms, color=palette[i], label=bin)
        bottoms += proportions
    ax.set_xticks(x, labels=df_wide.index, rotation=90, ha='right', rotation_mode='anchor')
    ax.set_xlabel("sample", labelpad=20)
    ax.set_ylabel(f"Proportion of {label}")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), title=f"{read_type} reads per cluster")
    fig.suptitle(f"Distribution of {label} by {read_type} reads per cluster")
    return fig


def main():
    # argument validation and processing
    args = parse_arguments()
    assert args.dpm_read or args.dpm_cluster or args.bpm_read or args.bpm_cluster, (
        "ERROR: At least one of the following arguments must be provided: --dpm_read, --dpm_cluster, "
        "--bpm_read, --bpm_cluster."
    )
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
        "antibody ID": 'category',
        "read type": 'category',
        "reads per cluster": int,
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
        df.groupby(["sample", "read type", "reads per cluster"], observed=True)
        ["number of clusters"].sum()
        .reset_index()
    )
    df_agg = pd.concat(
        (
            df_agg,
            (
                df_agg
                .groupby(["read type", "reads per cluster"], observed=True)
                ["number of clusters"].sum()
                .reset_index()
                .assign(sample="all")
            )
        ),
        axis=0,
        ignore_index=True
    ).sort_values(["sample", "read type", "reads per cluster"])
    df_agg['number of reads'] = df_agg['reads per cluster'] * df_agg['number of clusters']
    df_agg['reads per cluster, binned'] = pd.cut(
        df_agg["reads per cluster"],
        bins=[0, 0.5, 1, 5, 10, 20, 50, 100, 200, np.inf],
        labels=["0", "1", "2-5", "6-10", "11-20", "21-50", "51-100", "101-200", "201+"],
        include_lowest=True,
        right=True
    )
    mask_dpm = df_agg["read type"] == 'DPM'
    mask_bpm = df_agg["read type"] == 'BPM'

    if args.dpm_cluster:
        with PdfPages(args.dpm_cluster) as pdf:
            fig_ecdf = cluster_distribution_ecdf(df_agg.loc[mask_dpm], "DPM", "clusters", palette=palette)
            pdf.savefig(fig_ecdf, bbox_inches="tight")
            fig_stacked_bar = cluster_distribution_stacked_bar(df_agg.loc[mask_dpm], "DPM", "clusters", palette=palette)
            pdf.savefig(fig_stacked_bar, bbox_inches="tight")
    if args.dpm_read:
        with PdfPages(args.dpm_read) as pdf:
            fig_ecdf = cluster_distribution_ecdf(df_agg.loc[mask_dpm], "DPM", "reads", palette=palette)
            pdf.savefig(fig_ecdf, bbox_inches="tight")
            fig_stacked_bar = cluster_distribution_stacked_bar(df_agg.loc[mask_dpm], "DPM", "reads", palette=palette)
            pdf.savefig(fig_stacked_bar, bbox_inches="tight")
    if args.bpm_cluster:
        with PdfPages(args.bpm_cluster) as pdf:
            fig_ecdf = cluster_distribution_ecdf(df_agg.loc[mask_bpm], "BPM", "clusters", palette=palette)
            pdf.savefig(fig_ecdf, bbox_inches="tight")
            fig_stacked_bar = cluster_distribution_stacked_bar(df_agg.loc[mask_bpm], "BPM", "clusters", palette=palette)
            pdf.savefig(fig_stacked_bar, bbox_inches="tight")
    if args.bpm_read:
        with PdfPages(args.bpm_read) as pdf:
            fig_ecdf = cluster_distribution_ecdf(df_agg.loc[mask_bpm], "BPM", "reads", palette=palette)
            pdf.savefig(fig_ecdf, bbox_inches="tight")
            fig_stacked_bar = cluster_distribution_stacked_bar(df_agg.loc[mask_bpm], "BPM", "reads", palette=palette)
            pdf.savefig(fig_stacked_bar, bbox_inches="tight")


if __name__ == "__main__":
    main()
