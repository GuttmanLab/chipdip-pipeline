import argparse
import collections
import collections.abc
import contextlib
import itertools
import os
import sys
import typing

import numpy as np
import pandas as pd
import pysam
from tqdm.auto import tqdm

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers

"""
Label reads with antibody assignment via read group tags.
"""

UNASSIGNED_LABELS = ("filtered", "none", "uncertain", "ambiguous")


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Label reads with antibody assignment via read group tags."
    )
    parser.add_argument("config", metavar="config.txt", help="Barcode config file")
    parser.add_argument(
        "--config_type",
        choices=["splitcode", "bID"],
        required=True,
        help="Type of barcode config file: splitcode or bID (SPRITE).",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="in.bam",
        help=(
            "Path to input BAM file, where reads have been given read type (RT) and cluster barcode (CB) tags. "
            "If not provided, read from standard input."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="out.bam",
        help=(
            "Path to output BAM file, where DPM reads are given read group (RG) tags, and BPM reads are given "
            "read group (YG) and DPM count (YC) tags. If a hyphen (-), write to standard output. If not provided, "
            "no BAM output is written, and only cluster statistics are returned."
        )
    )
    parser.add_argument(
        "--output_reads_per_cluster",
        metavar="reads_per_cluster.tsv[.gz]",
        help=(
            "Histogram of reads per cluster, for each antibody and read type. 4-column tab-delimited table with "
            "header: antibody ID, read type (DPM or BPM), reads per cluster, number of clusters."
        )
    )
    parser.add_argument(
        "--output_bpm_max_rep",
        metavar="bpm_max_rep.tsv[.gz]",
        help=(
            "Histogram of counts or proportions for the maximally-represented BPM per cluster. 4-column tab-delimited"
            "table with header: antibody ID, metric, value, number of clusters, where metric can be 'count' or "
            "'proportion'."
        )
    )
    parser.add_argument(
        "--min_oligos",
        type=int,
        required=True,
        help="The minimum number of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--proportion",
        type=float,
        required=True,
        help="The minimum representation proportion of oligos needed to call a cluster",
    )
    parser.add_argument(
        "--max_size",
        type=int,
        required=True,
        help="The maximum cluster size to keep",
    )
    parser.add_argument(
        "--tag_barcode",
        default="CB",
        help="2-letter SAM/BAM tag name for the cluster barcode.",
    )
    parser.add_argument(
        "--tag_read_type",
        default="RT",
        help="2-letter SAM/BAM tag name for the read type.",
    )
    parser.add_argument(
        "--tag_read_group",
        default="RG",
        help="2-letter SAM/BAM tag name for the read group of DPM reads.",
    )
    parser.add_argument(
        "--tag_bpm_read_group",
        default="YG",
        help=(
            "2-letter SAM/BAM tag name for the read group of BPM reads. "
            "Set to the empty string to skip setting this tag."
        ),
    )
    parser.add_argument(
        "--tag_dpm_count",
        default="YC",
        help=(
            "2-letter SAM/BAM tag name for DPM counts per cluster for BPM reads. "
            "Set to the empty string to skip setting this tag."
        ),
    )
    parser.add_argument("-u", "--uncompressed", action="store_true", help="Uncompressed BAM output.")
    parser.add_argument("-t", "--threads", type=helpers.positive_int, default=1, help="Number of threads to use.")
    parser.add_argument("--progress", action="store_true", help="Display a progress bar.")
    args = parser.parse_args()
    return args


@contextlib.contextmanager
def open_bam(file: str | pysam.AlignmentFile | typing.IO | None, **kwargs) -> collections.abc.Generator[pysam.AlignmentFile | None]:
    """
    Context manager to open a BAM file if given as a string path, or yield the file object if already opened.

    Args:
    - file: Path to the BAM file or a pysam.AlignmentFile. If None, return None.
    - **kwargs: Additional arguments to pass to pysam.AlignmentFile.

    Yields: Opened BAM file or None.
    """
    if file is None:
        yield None
    elif isinstance(file, pysam.AlignmentFile):
        yield file
    else:
        with pysam.AlignmentFile(file, **kwargs) as bam_file:
            yield bam_file


def parse_splitcode_config(path_config: str) -> tuple[list[int], dict[str, str]]:
    """
    Parse a splitcode config file to get the following:
    1. number of tags in the barcode for each read orientation
    2. antibody ID names

    Args
    - path_config: Path to the splitcode config file

    Returns
    - num_tags: Number of tags in the barcode for each read orientation
    - antibody_IDs: Map from antibody ID names to antibody ID tag sequences; if multiple sequences correspond to the
      same antibody ID, then they are concatenated with a hyphen.
    """
    pos = 0
    with open(path_config, 'rt') as f:
        # skip blank lines, comments (#), and header options (@)
        line = f.readline()
        while line:
            if line.startswith(('#', '@')) or line.strip() == '':
                pos = f.tell()
            else:
                break
            line = f.readline()
        if line == '':
            raise ValueError(f"No valid lines found in barcode config file {path_config}")
        f.seek(pos)
        df = pd.read_csv(
            f,
            sep="\t",
            header=0,
            comment="#",
        ).rename(columns={'tag': 'tags', 'id': 'ids', 'group': 'groups', 'distance': 'distances', 'location': 'locations', 'sub': 'subs'})
        # From splitcode documentation: "the tags, ids, groups, distances, locations, and subs option can be specified either singular or plural"
        # --> map to plural column names for consistency
    if 'locations' in df.columns:
        # only look for antibody ID tag sequences in read 1
        mask_locations = df['locations'].astype(str).str.startswith(('0', '-1'))
    else:
        mask_locations = np.ones(len(df), dtype=bool)
    antibody_IDs = df.loc[df['ids'].str.startswith("BEAD_") & mask_locations].groupby('ids')['tags'].agg(lambda x: '-'.join(x)).to_dict()
    antibody_IDs = {k[5:]: v for k, v in antibody_IDs.items()} # remove "BEAD_" prefix
    num_tags = (
        df.loc[df['exclude'].astype(str) == '-']
        .assign(read=lambda tb: tb['locations'].str[0])
        .groupby('read')['groups']
        .unique().map(len)
        .sort_index()
        .tolist()
    )
    return num_tags, antibody_IDs


def parse_bID_config(path_config: str) -> tuple[int, dict[str, str]]:
    """
    Parse a SPRITE barcode config file to get the following:
    1. number of tags in the barcode, excluding DPM and BEAD_ tags
    2. antibody ID names

    Args
    - path_config: Path to the barcode config file

    Returns
    - num_tags: Number of tags in the barcode, excluding DPM and BEAD_ tags.
    - antibody_IDs: Map from antibody ID names to antibody ID tag sequences; if multiple sequences correspond to the
      same antibody ID, then they are concatenated with a hyphen.
    """
    num_tags = 0
    read1_seen = False
    read2_seen = False
    antibody_IDs = dict()
    with open(path_config, 'rt') as f:
        for line in f:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            if line[:5].upper() == 'READ1':
                if read1_seen:
                    raise ValueError("Multiple READ1 lines found in config file.")
                num_tags += line.count("ODD") + line.count("EVEN") + line.count("Y")
                read1_seen = True
            elif line[:5].upper() == 'READ2':
                if read2_seen:
                    raise ValueError("Multiple READ2 lines found in config file.")
                num_tags += line.count("ODD") + line.count("EVEN") + line.count("Y")
                read2_seen = True
            else:
                entries = line.split('\t')
                if entries[0].upper() == 'DPM':
                    assert len(entries) == 4, \
                        "DPM lines should have 4 columns: 'DPM', tag name, tag sequence, tag error tolerance"
                    if entries[1].startswith("BEAD_"):
                        name = entries[1][5:] # remove "BEAD_" prefix
                        seq = entries[2]
                        if name in antibody_IDs:
                            seqs = antibody_IDs[name].split('-')
                            if seq not in seqs:
                                antibody_IDs[name] += '-' + seq
                        else:
                            antibody_IDs[name] = seq
    return num_tags, antibody_IDs


def construct_read_group_header(
    file_in: str | pysam.AlignmentFile,
    rgids_to_bcs: collections.abc.Mapping[str, str | None]
) -> dict:
    """
    Add read group lines to a BAM file header.

    Args
    - file_in: Path to input BAM file, or a pysam.AlignmentFile object.
    - rgids_to_bcs: map from read group IDs to read group barcode sequences.
        Read group barcode sequences can be None or the empty string if not available.
        If there are several barcodes for each read group, the barcodes should be concatenated with hyphens ('-').

    Returns
    - bam_header: Updated BAM header
    """
    with open_bam(file_in) as bam:
        bam_header = bam.header.to_dict()
        read_groups = []
        for rgid, bc in rgids_to_bcs.items():
            if bc is None or bc == '':
                read_groups.append({"ID": rgid})
            else:
                read_groups.append({"ID": rgid, "BC": bc})
        bam_header["RG"] = read_groups
        return bam_header


def assign_label(
    cluster: collections.abc.Collection[pysam.AlignedSegment],
    thresh_min_oligos: int,
    thresh_proportion: float,
    max_size: int | None = None,
    tag_read_type: str = "RT",
) -> tuple[str, int, int, int, float]:
    """
    Assign label to a cluster of reads.

    Args
    - cluster: Collection of reads, assumed to be deduplicated.
    - min_oligos: The minimum number of oligos needed to call a cluster.
    - proportion: Minimum representation proportion of oligos needed to call a cluster.
    - max_size: Maximum number of DNA (DPM) reads allowed per cluster.
        If the maximum size is exceeded, then the cluster is assigned the label "filtered".
        If None, no maximum size is enforced.
    - tag_read_type: 2-letter code to use for the read type tag.

    Returns
    - label: Antibody ID or an unassigned label ("filtered", "none", "uncertain", or "ambiguous").
    - n_DPM_reads: Number of DNA (DPM) reads in the cluster.
    - n_skipped: Number of reads skipped due to unexpected read type.
    - count: Number of BPM reads from the maximally represented antibody ID in the cluster.
    - proportion: Number of BPM reads from the maximally represented antibody ID in the cluster.
    """
    n_DPM_reads = 0
    n_skipped = 0
    antibody_IDs_counter = collections.Counter()
    for read in cluster:
        read_type = read.get_tag(tag_read_type)
        if not isinstance(read_type, str):
            print(f"Unexpected read type: {read_type}", file=sys.stderr)
            n_skipped += 1
            continue
        if read_type.startswith("DPM"):
            n_DPM_reads += 1
        elif read_type.startswith("BEAD_"):
            antibody_IDs_counter[read_type[5:]] += 1
        else:
            print(f"Unexpected read type: {read_type}", file=sys.stderr)
            n_skipped += 1
    total = antibody_IDs_counter.total()
    if total == 0:
        label = "none"
        count = 0
        proportion = 0.0
    else:
        candidate, count = antibody_IDs_counter.most_common()[0]
        proportion = count / total
        if max_size is not None and n_DPM_reads > max_size:
            label = "filtered"
        else:
            if count < thresh_min_oligos:
                label = "uncertain"
            elif proportion < thresh_proportion:
                label = "ambiguous"
            else:
                label = candidate
    return label, n_DPM_reads, n_skipped, count, proportion


def label_bam_file(
    file_in: str | typing.IO,
    file_out: str | typing.IO | None,
    antibody_IDs: collections.abc.Mapping[str, str | None],
    min_oligos: int,
    proportion: float,
    max_size: int | None = None,
    tag_barcode: str = "CB",
    tag_read_type: str = "RT",
    tag_read_group: str = "RG",
    tag_bpm_read_group: str | None = "YG",
    tag_dpm_count: str | None = "YC",
    uncompressed: bool = False,
    threads: int = 1,
    progress: bool = False
) -> tuple[
    collections.Counter[tuple[str, str, int]],
    collections.Counter[tuple[str, str, float]],
    int
]:
    """
    Compute cluster statistics, and add antibody label to individual reads in a BAM file as read group tags.

    Args
    - file_in: Path or file object to input BAM file
    - file_out: Path or file object to write labeled BAM file. If None, no output is written, and only cluster
        statistics are returned.
    - antibody_IDs: Map from antibody ID names to antibody ID tag sequences; if multiple sequences correspond to the
      same antibody ID, then they should be concatenated with a hyphen. Tag sequences can be None or the empty string if
      not available.
    - min_oligos: The minimum number of oligos needed to call a cluster.
    - proportion: Minimum representation proportion of oligos needed to call a cluster.
    - max_size: Maximum number of DNA (DPM) reads allowed per cluster.
        If the maximum size is exceeded, then the cluster is assigned the label "filtered".
        If None, no maximum size is enforced.
    - tag_barcode: 2-letter SAM/BAM tag name for the cluster barcode.
    - tag_read_type: 2-letter SAM/BAM tag name for the read type.
    - tag_read_group: 2-letter SAM/BAM tag name for the read group of DPM reads.
    - tag_bpm_read_group: 2-letter SAM/BAM tag name for the read group of BPM reads.
        Set to the empty string to skip setting this tag.
    - tag_dpm_count: 2-letter SAM/BAM tag name for DPM counts per cluster for BPM reads.
        Set to the empty string to skip setting this tag.
    - uncompressed: Write uncompressed BAM output.
    - threads: Number of threads to use for writing the output BAM. Only applicable if uncompressed is False.
    - progress: Display a progress bar.

    Returns
    - reads_per_cluster: Map from (antibody ID, read type, reads per cluster) to number of clusters
    - bpm_max_rep: Map from (antibody ID, metric, value) to number of clusters, where metric is 'count' or 'proportion',
        and value is the count or proportion of BPM reads from the maximally represented antibody ID in the cluster,
        with value rounded to 2 decimal places for proportions.
    - skipped: Number of reads skipped due to unexpected read type, unexpected antibody IDs, or other errors.
    """
    antibody_IDs_and_unassigned = set(antibody_IDs.keys()) | set(UNASSIGNED_LABELS)
    mode_out = "wb"
    if uncompressed:
        mode_out = "wbu"
        threads = 1

    # cluster statistics counters
    # - reads_per_cluster: map from (antibody ID, read type, reads per cluster) to number of clusters
    # - bpm_max_rep: map from (antibody ID, metric, value) to number of clusters, where metric is 'count' or
    #     'proportion', and value is the count or proportion of BPM reads from the maximally represented antibody ID in
    #     the cluster, with value rounded to 2 decimal places for proportions.
    reads_per_cluster = collections.Counter()
    bpm_max_rep = collections.Counter()

    skipped = 0
    with pysam.AlignmentFile(file_in, mode="r") as bam_in:
        header = construct_read_group_header(bam_in, dict(antibody_IDs) | {label: None for label in UNASSIGNED_LABELS})
        with open_bam(file_out, mode=mode_out, header=header, threads=threads) as bam_out:
            for barcode, reads in itertools.groupby(
                tqdm(bam_in.fetch(until_eof=True), disable=not progress),
                key=lambda read: read.get_tag(tag_barcode)
            ):
                reads = list(reads)
                antibody_ID, n_DPM_reads, n_other, max_rep_count, max_rep_proportion = assign_label(
                    reads,
                    thresh_min_oligos=min_oligos,
                    thresh_proportion=proportion,
                    max_size=max_size,
                    tag_read_type=tag_read_type
                )

                # cluster statistics
                reads_per_cluster[(antibody_ID, "total", len(reads) - n_other)] += 1
                reads_per_cluster[(antibody_ID, "DPM", n_DPM_reads)] += 1
                reads_per_cluster[(antibody_ID, "BPM", len(reads) - n_DPM_reads - n_other)] += 1
                bpm_max_rep[(antibody_ID, "count", max_rep_count)] += 1
                bpm_max_rep[(antibody_ID, "proportion", round(max_rep_proportion, 2))] += 1

                if antibody_ID not in antibody_IDs_and_unassigned:
                    print(
                        f"Warning: Unexpected antibody ID '{antibody_ID}' found in cluster '{barcode}'",
                        file=sys.stderr
                    )
                    skipped += len(reads)
                    continue

                if bam_out is not None:
                    for read in reads:
                        read_type = read.get_tag(tag_read_type)
                        if not isinstance(read_type, str):
                            print(f"Unexpected read type ({read_type}) in read {read.query_name}.", file=sys.stderr)
                            skipped += 1
                            continue
                        if read_type.startswith("DPM"):
                            read.set_tag(tag_read_group, antibody_ID, replace=True)
                        elif read_type.startswith("BEAD_"):
                            # oligo reads are assigned YG and YC tags for antibody ID and count of DPM reads in the
                            # cluster; the RG tag is removed. Using different antibody ID tags for DPM (RG) and BPM (YG)
                            # reads allows subsequent samtools split to separately output DPM and BPM reads into their
                            # own files.
                            read.set_tag(tag_read_group, None, replace=True)
                            if tag_bpm_read_group:
                                read.set_tag(tag_bpm_read_group, antibody_ID, replace=True)
                            if tag_dpm_count:
                                read.set_tag(tag_dpm_count, n_DPM_reads, replace=True)
                        else:
                            print(f"Unexpected read type ({read_type}) in read {read.query_name}.", file=sys.stderr)
                            skipped += 1
                            continue
                        bam_out.write(read)

    # check that aggregating counts in reads_per_cluster matches the counts in splitbam_counts
    # for (antibody_ID, read_type, reads_per_cluster), count in reads_per_cluster.items():
    #     if read_type == "DPM":
    #         assert splitbam_counts[antibody_ID] >= count, (
    #             f"Mismatch in counts for {antibody_ID} DPM reads: {splitbam_counts[antibody_ID]} < {count

    # print("Reads written:", splitbam_counts.total(), file=sys.stderr)
    print("Reads with an error not written out:", skipped, file=sys.stderr)
    return reads_per_cluster, bpm_max_rep, skipped


def main():
    args = parse_args()

    # check that at least 1 output is specified
    if args.output is None and args.output_reads_per_cluster is None and args.output_bpm_max_rep is None:
        raise ValueError(
            "At least one of --output, --output_reads_per_cluster, or --output_bpm_max_rep must be specified."
        )

    file_out = None
    if args.output is not None:
        file_out = args.output if args.output != '-' else sys.stdout

    if args.config_type == "splitcode":
        _, antibody_IDs = parse_splitcode_config(args.config)
    else:  # args.config_type == "bID"
        _, antibody_IDs = parse_bID_config(args.config)

    reads_per_cluster, bpm_max_rep, skipped = label_bam_file(
        args.input if args.input else sys.stdin,
        file_out,
        antibody_IDs,
        args.min_oligos,
        args.proportion,
        args.max_size,
        tag_barcode=args.tag_barcode,
        tag_read_type=args.tag_read_type,
        tag_read_group=args.tag_read_group,
        tag_bpm_read_group=args.tag_bpm_read_group,
        tag_dpm_count=args.tag_dpm_count,
        uncompressed=args.uncompressed,
        threads=args.threads,
        progress=args.progress,
    )

    if args.output_reads_per_cluster:
        df_reads_per_cluster = (
            pd.Series(reads_per_cluster, name='number of clusters')
            .rename_axis(['antibody ID', 'read type', 'reads per cluster'])
            .reset_index()
            .sort_values(['antibody ID', 'read type', 'reads per cluster'])
        )
        df_reads_per_cluster.to_csv(args.output_reads_per_cluster, sep='\t', index=False, header=True)
        # to convert to splitbam counts:
        # (
        #     df_reads_per_cluster
        #     .loc[df_reads_per_cluster['read type'] == "DPM"]
        #     .assign(**{'number of reads': lambda tb: tb['reads per cluster'] * tb['number of clusters']})
        #     .groupby('antibody ID')['number of reads'].sum()
        # )
    if args.output_bpm_max_rep:
        (
            pd.Series(bpm_max_rep, name='number of clusters')
            .rename_axis(['antibody ID', 'metric', 'value'])
            .reset_index()
            .sort_values(['antibody ID', 'metric', 'value'])
            .to_csv(args.output_bpm_max_rep, sep='\t', index=False, header=True)
        )


if __name__ == "__main__":
    main()
