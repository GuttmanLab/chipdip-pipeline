"""
Rename and select chromosomes in a FASTA or BAM file.
"""

import argparse
import copy
import gzip
import io
import os
import shutil
import subprocess
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers
import pysam


def main():
    """
    Parse arguments and execute the following behavior based on the arguments as follows:

    chrom_map | try_symlink | output | behavior
    --------- | ----------- | ------ | --------
    None      | True, False | None   | Write input to standard out
    None      | True        | <path> | Create symbolic link from input to output <path>;
                                       output <path> must not point to an existing file
    None      | False       | <path> | Copy input to output <path>
    <path>    | True, False | None   | Rename/filter chromosomes, write to standard out
    <path>    | True, False | <path> | Rename/filter chromosomes, write to <path>
    """
    args = parse_arguments()

    # support standard input
    if args.input == "-":
        args.input = sys.stdin.buffer
        args.try_symlink = False

    if args.chrom_map is None:
        if args.output is None:
            # write input to standard out
            if args.input is sys.stdin.buffer:
                shutil.copyfileobj(sys.stdin.buffer, sys.stdout.buffer)
            else:
                with open(args.input, "rb") as f:
                    shutil.copyfileobj(f, sys.stdout.buffer)
        else:
            # try symlink if requested
            if args.try_symlink:
                try:
                    os.symlink(os.path.abspath(args.input), args.output)
                    return
                except Exception as err:
                    if not args.quiet:
                        print(
                            f"Error upon attempt to create a symbolic link from {args.output} to {args.input}:", err
                        )
            # if symlink fails or symlink not requested, copy input to output
            if args.input is sys.stdin.buffer:
                with open(args.output, "wb") as f:
                    shutil.copyfileobj(sys.stdin.buffer, f)
            else:
                shutil.copyfile(args.input, args.output)
    else:
        chrom_map = helpers.parse_chrom_map(args.chrom_map)
        if args.fasta:
            filter_fasta(
                args.input,
                args.output,
                chrom_map
            )
        else:
            filter_reads(
                args.input,
                args.output,
                chrom_map,
                try_symlink=args.try_symlink,
                sort=args.sort,
                no_PG=args.no_PG,
                threads=args.threads,
                verbose=not args.quiet,
            )



def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Rename and select chromosomes in a FASTA or BAM file. "
            "For BAM files, keep only reads aligned to selected chromosomes, and "
            "reorder the chromosomes in the BAM file header."
        )
    )
    parser.add_argument(
        "input",
        metavar="in.bam|in.fasta(.gz)|-",
        help=(
            "Input file. Assumed to be BAM format unless the -f/--fasta flag is used. "
            "Use '-' for standard input."
        )
    )
    parser.add_argument(
        "-f", "--fasta",
        action="store_true",
        help="Input file is FASTA format."
    )
    parser.add_argument(
        "-o", "--output",
        metavar="out.bam|out.fasta",
        help="Output file. If not provided, writes to standard out."
    )
    parser.add_argument(
        "-c", "--chrom_map",
        metavar="PATH",
        help="Chromosome name map file."
    )
    parser.add_argument(
        "--try-symlink",
        action="store_true",
        help=(
            "If -c/--chrom_map is not specified (i.e., no changes are required), "
            "try to make a symbolic link from input to output."
        )
    )
    parser.add_argument(
        "-t", "--threads",
        type=helpers.positive_int,
        default=1,
        metavar="#",
        help="Number of threads to use for compressing/decompressing BAM files."
    )
    parser.add_argument(
        "--sort",
        type=str,
        choices=("auto", "true", "false"),
        default="auto",
        help=(
            "(Only relevant if -c/--chrom_map is specified and the input is a BAM file) "
            "Whether to sort the processed BAM file. "
            "If 'auto', will only re-sort if the order of chromosomes given in the chromosome name map "
            "is different than the existing chromosome order."
        )
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Do not print the number of discarded reads and output reads to standard error."
    )
    parser.add_argument(
        "--no-PG",
        action="store_true",
        help="Do not add @PG line to the header of the output file if sorting BAM file."
    )
    return parser.parse_args()


def reheader(old_header, chrom_map):
    """
    Rename and reorder reference sequence names. Update "@HD SO:" and "@HD SS:" header tags
    as appropriate: if reference sequences are reordered, then set "@HD SO:unsorted" and
    remove the "@HD SS:" tag.

    Args
    - old_header: dict (str -> str)
        Representation of SAM/BAM file header of existing BAM file as a Python
        dictionary, as used by pysam. For example, generated by
        `pysam.AlignmentFile.header.to_dict()`
    - chrom_map: dict (str -> str)
        Map from old reference sequence names to new reference sequence names.
        The order of entries in this dictionary defines the new alignment sorting order.

    Returns
    - new_header: dict
        Representation of SAM/BAM file header as a Python dictionary, as used by pysam
    - old_to_new_refID: dict (int -> int)
        Map from old reference ID to new reference ID.
        Note that BAM reference IDs are 0-indexed.
    - retains_sorting: bool
        Whether the order of entries in chrom_map is consistent with existing order
        of chromosomes in the old header.

    Note: The current implementation of this function takes advantage of dict features
    defined/introduced in Python versions >= 3.9.
    """
    # copy existing BAM header to a new header dict without reference sequences
    new_header = copy.deepcopy(old_header)
    old_name_to_entries = {entry["SN"]: (refID, entry) for (refID, entry) in enumerate(new_header["SQ"])}
    new_header["SQ"] = []

    old_to_new_refID = {}
    for new_refID, (old_name, new_name) in enumerate(chrom_map.items()):
        old_refID, entry = old_name_to_entries[old_name]
        entry["SN"] = new_name
        new_header["SQ"].append(entry)
        old_to_new_refID[old_refID] = new_refID

    # check if the order of the chromosomes have changed
    old_refID_order = list(old_to_new_refID.keys())
    if len(chrom_map) == 1:
        retains_sorting = True
    else:
        retains_sorting = all(a < b for a, b in zip(old_refID_order, old_refID_order[1:]))

    if not retains_sorting:
        new_header["HD"]["SO"] = "unsorted"
        new_header["HD"].pop("SS", None)
    return new_header, old_to_new_refID, retains_sorting


def filter_reads(
    path_bam_in, path_bam_out, chrom_map, try_symlink=False, sort="auto", no_PG=False, threads=1, verbose=True
):
    """
    Discard reads that do not map to the specified chromosomes, and generate a new header for only the specified
    chromosomes.

    Args
    - path_bam_in: str, io.TextIOWrapper, io.BufferedReader
        Path to input BAM file or sys.stdin or sys.stdin.buffer
    - path_bam_out: str or None
        Path to output BAM file. If None, writes to standard out.
    - chrom_map: dict (str -> str)
        Map from old reference sequence names to new reference sequence names.
        The order of entries in this dictionary defines the new alignment sorting order.
    - try_symlink: bool. default=False
        If no renaming, filtering, or sorting is necessary, try to create a symbolic link from
        path_bam_in to path_bam_out.
    - sort: ('auto', 'true', or 'false'). default='auto'
        Whether to coordinate sort the reads before writing to path_bam_out.
    - no_PG: bool. default=False
        Whether to add @PG line to to the header of the output file if sorting BAM file.
    - threads: int. default=1
        Number of threads to use for compressing/decompressing BAM files
    - verbose: bool. default=True
        Print the number of discarded reads and the number of output reads.
    """
    count_discard = 0
    count_out = 0
    with pysam.AlignmentFile(path_bam_in, "rb", threads=threads) as file_bam_in:
        old_header = file_bam_in.header.to_dict()
        new_header, old_to_new_refID, retains_sorting = reheader(old_header, chrom_map)

        if (old_header == new_header) and (sort != "true") and (path_bam_out is not None) and try_symlink:
            assert retains_sorting
            try:
                os.symlink(os.path.abspath(path_bam_in), path_bam_out)
                return
            except Exception as err:
                if verbose:
                    print(f"Error upon attempt to create a symbolic link from {path_bam_in} to {path_bam_out}:", err)

        def process_reads(output_stream):
            nonlocal count_discard, count_out
            for read in file_bam_in.fetch(until_eof=True):
                if read.reference_name in chrom_map:
                    read.reference_id = old_to_new_refID[read.reference_id]
                    # if reference sequence name of paired read is not in chrom_map, set RNEXT
                    # to be "*"
                    read.next_reference_id = old_to_new_refID.get(read.next_reference_id, -1)
                    output_stream.write(read)
                    count_out += 1
                else:
                    count_discard += 1

        if sort == "true" or ((sort == "auto") and retains_sorting is False):
            sort_cmd = ["samtools", "sort"]
            if path_bam_out is not None:
                sort_cmd.extend(["-o", path_bam_out])
                stdout = None
            else:
                stdout = sys.stdout.buffer
            if no_PG:
                sort_cmd.append("--no-PG")
            with subprocess.Popen(sort_cmd, stdin=subprocess.PIPE, stdout=stdout) as popen_samtools:
                with pysam.AlignmentFile(
                    popen_samtools.stdin,
                    "wb",
                    header=new_header,
                    threads=threads
                ) as file_bam_out:
                    process_reads(file_bam_out)
            if verbose:
                print(popen_samtools, file=sys.stderr)
        else:
            with pysam.AlignmentFile(
                path_bam_out if path_bam_out is not None else sys.stdout.buffer,
                "wb",
                header=new_header,
                threads=threads,
            ) as file_bam_out:
                process_reads(file_bam_out)
        sys.stdout.flush()

    if verbose:
        print("Discarded reads:", count_discard, file=sys.stderr)
        print("Written out reads:", count_out, file=sys.stderr)


def filter_fasta(path_in, path_out, chrom_map) -> None:
    """
    Rename and select chromosomes from a FASTA file according to a chromosome name map.

    Args
    - path_in: str or file object
        Path to input FASTA file or file object.
        Support normal text or gzip-compressed text.
    - path_out: str, file object, or None
        Path to output FASTA file or file object. If None, writes to standard out.
    - chrom_map: dict (str -> str)
        Map from old reference sequence names to new reference sequence names.
    """
    # convert path_in to a text IO stream if necessary
    if isinstance(path_in, io.IOBase) and not isinstance(path_in, io.TextIOBase):
        first_two_bytes = path_in.peek(2)[:2]
        if first_two_bytes == helpers.GZIP_MAGIC_NUMBER:
            f = gzip.GzipFile(fileobj=path_in, mode="rb")
        else:
            f = path_in
        return filter_fasta(io.TextIOWrapper(f, encoding="utf-8"), path_out, chrom_map)
    elif not isinstance(path_in, io.IOBase):
        with helpers.file_open(path_in, mode="rt") as f:
            return filter_fasta(f, path_out, chrom_map)

    # convert path_out to a text IO stream if necessary
    if path_out is None:
        path_out = sys.stdout
    if not isinstance(path_out, io.IOBase):
        if path_out.endswith('.gz'):
            with gzip.open(path_out, "wt") as f:
                return filter_fasta(path_in, f, chrom_map)
        else:
            with open(path_out, "wt") as f:
                return filter_fasta(path_in, f, chrom_map)

    include = False
    for line in path_in:
        if line.startswith('>'):
            old_name = helpers.REGEX_RNAME.search(line[1:]).group()
            if old_name in chrom_map:
                path_out.write('>' + chrom_map[old_name] + '\n')
                include = True
            else:
                include = False
        elif include:
            path_out.write(line)

if __name__ == "__main__":
    main()
