"""
Demultiplex a BAM file into multiple BAM files based on sample information encoded in the barcode.

Notes
- The word "tag" is overloaded: infer from context whether it refers to a SAM/BAM tag (e.g., CB) or a tag within a
  barcode (e.g., "O1" for an index-1 Odd tag in a SPRITE cluster barcode).

The config file is inspired by splitcode's config file and is formatted as follows:
- Lines starting with # are comments and ignored.
- Lines that are empty or contain only whitespace are ignored.
- The file consists of 2 parts, always in this order:
  - Rules
    - Lines starting with @keep define selectors (tag combinations) corresponding to each sample.
      - @keep lines are formatted as: @keep <tab> selector1,selector2,... <tab> sample_name
        - sample_name can contain whitespace but not tab characters.
    - Lines starting with @remove define selectors to be discarded.
      - @remove lines are formatted as: @remove <tab> selector1,selector2,...
    - selectors are comma-delimited lists of tag IDs or group names enclosed in curly braces: {id}|{{group}}[,{id}|{{group}},...]
      - If a tag ID is used, the read must match that tag.
      - If a group name is used, the read must match at least one tag in that group.
      - A read is assigned to a sample if it matches at least one tag from each selector in the @keep line for that sample.
  - Tag definitions
    - The first line of this section is a header with columns "id", "tag", "location", and optionally "group".
    - Each subsequent line defines a tag definition:
      - id: unique identifier for the tag definition (each line must have a unique id value; should not start with '@'
      - tag: the tag name in a barcode string (e.g., "O1" for Odd tag 1)
      - location: the SAM/BAM tag where the barcode is stored, along with the 0-based index of the tag within the barcode string,
        formatted as "FIELD:INDEX" (e.g., "CB:2" for the 3rd position in the CB tag)
      - group: (optional) group name to which this tag belongs; can be referenced in selectors using {{group}}. Use a '-'
        character to indicate no group.

Example config file:

@keep	{tag1},{{group1}}	sample1
@keep	{{group2}}	sample2
@remove	{tag1},{{group2}}

id	tag	location	group
tag1	O1	CB:0	group0
tag2	O2	CB:0	group0
tag3	DPM1	RT:1	group1
tag4	DPM2	RT:1	group1
tag5	DPM3	RT:1	group2
tag6	DPM4	RT:1	group2
"""
import argparse
import collections
import contextlib
import csv
import re
import pdb
import sys
import typing

import pysam
from tqdm.auto import tqdm

COLS_HEADER_REQUIRED = {"id", "tag", "location"}
COLS_HEADER_ALL = COLS_HEADER_REQUIRED | {"group"}

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Demultiplex a BAM file into multiple BAM files based on sample information encoded in the barcode."
    )
    parser.add_argument(
        '-i', '--input',
        help='Input BAM file. If not provided, read from standard input.'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output BAM file path. At least one of --output or --output_pattern must be provided.'
    )
    parser.add_argument(
        '--output_pattern',
        help=(
            'Format string for output BAM paths when demultiplexing into separate BAM files per sample; must include '
            'the "{sample}" field. BAM files will be created even for samples with no matching reads. At least one of '
            "--output or --output_pattern must be provided. Example: 'samples/{sample}.bam'"
        )
    )
    parser.add_argument(
        '-u', '--unassigned',
        help=(
            'Output BAM file path for reads not matching any sample and not discarded by @remove directives. If not '
            'provided, unmatched reads are discarded.'
        )
    )
    parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to TSV file defining tags, groups, and rules for demultiplexing.'
    )
    parser.add_argument(
        '--allow_multiple',
        action='store_true',
        help='Allow reads to match multiple samples; if not set, only the first matching sample will be used.'
    )
    parser.add_argument(
        '--sample_tag',
        help=(
            '2-letter SAM/BAM tag to use for storing sample information; will overwrite existing values. If "RG", '
            'existing @RG header lines in the SAM/BAM file will be discarded, and new header lines will be added. If '
            'allow_multiple is set, matching sample names will be semicolon-delimited. Unassigned reads do not get a '
            'sample tag. If not provided, no tag will be added to reads.'
        )
    )
    parser.add_argument(
        '--delimiter',
        default='.',
        help='Delimiter to split the barcode string. If set to the empty string "", match the entire barcode string.'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for output BAM compression.'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Print status messages and show progress bar.'
    )
    return parser.parse_args()


class TagDef:
    """
    Corresponds to a single line (tag definition) in the sample config file.
    """
    def __init__(self, tag_id: str, tag: str, location: str, group: str | None = None):
        self.id = tag_id
        self.tag = tag
        self.group = group if group != '-' else None
        self.field, idx = location.split(":")
        self.index = int(idx)

    def matches(self, read: pysam.AlignedSegment, delimiter: str = ".") -> bool:
        """
        Check if a read matches this tag definition.

        Args
        - read: read
        - delimiter: delimiter to split the barcode string; if the empty string, match the entire string

        Returns: True if the read matches the tag definition, False otherwise.
        """
        try:
            val = read.get_tag(self.field)
        except KeyError:
            return False

        if delimiter:
            parts = val.split(delimiter)
            return self.index < len(parts) and parts[self.index] == self.tag
        else:
            return val == self.tag


class Rule:
    """
    Corresponds to a single keep or remove rule in the sample config file.
    """
    def __init__(self, selectors: list[set[str]], sample=None):
        self.selectors = selectors
        self.sample = sample

    def matches(self, present: set[str]) -> bool:
        """
        Given a set of tag ids present in a read, check if the read matches this rule.
        Example usage:
        - self.selectors = [{"tag1"}, {"tag2", "tag3"}]
        - present = {"tag1", "tag3"}
        - returns True because the read has at least one tag from each selector set.
        """
        return all(any(t in present for t in sel) for sel in self.selectors)


def parse_selectors(expr: str, groups: dict[str, set[str]]) -> list[set[str]]:
    """
    Parse a selector expression into a list of selectors. A selector is a set of tag ids.

    Args
    - expr: comma-delimited list of selectors, where each selector is {id} or {{group}}.
        Example: "{tag1},{{group1}}"
    - groups: map of group name -> set of tag ids

    Returns: list of sets of tag ids corresponding to each selector.

    Example: Consider a line in the config file:
        @keep	{tag1},{{group1}}	sample1
    If group1 contains tags tag2 and tag3, then:
        arguments:
          - expr = "{tag1},{{group1}}"
          - groups = {"group1": {"tag2", "tag3"}}
        output:
          - [{"tag1"}, {"tag2", "tag3"}]
    """
    selectors = []
    for token in expr.split(","):
        selector = token.strip()
        if selector.startswith("{{") and selector.endswith("}}"):
            selectors.append(groups[selector[2:-2]])
        elif selector.startswith("{") and selector.endswith("}"):
            selectors.append({selector[1:-1]})
        else:
            raise ValueError(f"Invalid selector: {selector}")
    return selectors


def parse_config(path: str) -> tuple[dict[str, TagDef], dict[str, TagDef], list[Rule], list[Rule]]:
    """
    Parse the sample config file.

    Returns
    - tag_defs: dict mapping tag id -> TagDef; only includes tag definitions used in any rule
    - unused_tag_defs: dict mapping tag id -> TagDef; only includes tag definitions not used in any rule
    - keep_rules: list of Rule objects for @keep directives
    - remove_rules: list of Rule objects for @remove directives
    """
    tag_defs = dict()           # dict[str, TagDef]: map tag id -> TagDef
    groups = collections.defaultdict(set)   # dict[str, set[str]]: map group name -> set of tag ids
    keep_rules = []             # list[Rule]
    remove_rules = []           # list[Rule]

    with open(path) as f:
        lines = [l.strip() for l in f if (s := l.strip()) and s[0] != "#"] # walrus operator := requires Python 3.8+

    # parse tag definitions
    reader = csv.DictReader((l for l in lines if not l.startswith("@")), delimiter="\t")
    for row in reader:
        td = TagDef(row["id"], row["tag"], row["location"], row.get("group"))
        if td.id in tag_defs:
            raise ValueError(f"Duplicate tag id in config file: {td.id}")
        tag_defs[td.id] = td
        if td.group and td.group != '-':
            groups[td.group].add(td.id)

    # parse keep/remove rules
    for line in lines:
        if line.startswith("@keep"):
            expr, sample = line.split("\t")[1:3]
            selectors = parse_selectors(expr, groups)
            keep_rules.append(Rule(selectors, sample))
        elif line.startswith("@remove"):
            expr = line.split("\t", maxsplit=1)[1]
            selectors = parse_selectors(expr, groups)
            remove_rules.append(Rule(selectors))
        else:
            continue

    defined_tag_ids = set(tag_defs.keys()) # all defined tag ids
    used_tag_ids = set() # tag ids used in any rule
    for rule in keep_rules + remove_rules:
        for selector in rule.selectors:
            used_tag_ids.update(selector)

    # validate that all tag ids used in rules are defined
    for tag_id in used_tag_ids - defined_tag_ids:
        raise ValueError(f"Undefined tag id in config file: {tag_id}")

    # move unused tag definitions to unused_tag_defs
    unused_tag_defs = {tag_id: tag_defs.pop(tag_id) for tag_id in defined_tag_ids - used_tag_ids}

    return tag_defs, unused_tag_defs, keep_rules, remove_rules


def new_header_with_rg(original_header: dict[str, dict | list[dict]], samples: set[str]) -> dict:
    """
    Create a new SAM/BAM header with @RG lines for each sample.

    Args
    - original_header: original SAM/BAM header, as returned by pysam.AlignmentFile.header.to_dict()
        two-level dictionary where the first level contains the record (HD, SQ, etc) and the second level contains the
        fields (VN, LN, etc).
    - samples: set of sample names

    Returns: a copy of the original header with updated @RG lines
    """
    new_header = original_header.copy()
    if len(samples) > 0:
        new_header['RG'] = [{'ID': sample} for sample in samples]
    return new_header


def demultiplex_bam(
    config: str,
    bam_in: str | typing.IO,
    bam_out: str | typing.IO | None = None,
    bam_unassigned: str | typing.IO | None = None,
    output_pattern: str | None = None,
    allow_multiple: bool = True,
    sample_tag: str | None = None,
    delimiter: str = ".",
    threads: int = 1,
    verbose: bool = False,
) -> dict[str, int]:
    """
    Demultiplex a BAM file based on sample information encoded in barcodes.

    Args
    - config: Path to TSV file defining tags, groups, and rules for demultiplexing.
    - bam_in: Path to input BAM file or file-like object (e.g., sys.stdin).
    - bam_out: Path to output BAM file or file-like object (e.g., sys.stdout) for sample-assigned reads. If None,
        assigned reads are not written to a single output file.
    - bam_unassigned: Path to output BAM file for unassigned reads. If None, unassigned reads are discarded.
    - output_pattern: Format string for output BAM paths when demultiplexing into separate BAM files per sample; must
        include the "{sample}" field. BAM files will be created even for samples with no matching reads.
    - allow_multiple: Allow reads to match multiple samples; if False, only the first matching sample will be used.
    - sample_tag: 2-letter SAM/BAM tag to use for storing sample information; will overwrite existing values. If "RG",
        existing @RG header lines in the SAM/BAM file will be discarded, and new header lines will be added. If allow_multiple
        is set, matching sample names will be semicolon-delimited. If None, no tag will be added to reads.
    - delimiter: Delimiter to split the barcode string; if the empty, match the entire string.
    - threads: Number of threads to use for writing the output BAM files.
    - verbose: Print status messages and show progress bar.

    Returns: dictionary of counts
    - keep: Number of sample-assigned reads matching @keep rules
    - remove: Number of reads removed by @remove rules
    - unassigned: Number of reads not matching any sample and not removed by @remove rules
    """
    if verbose:
        print("Parsing config file...", file=sys.stderr)
    tag_defs, _, keep_rules, remove_rules = parse_config(config)
    samples = {r.sample for r in keep_rules}
    if verbose:
        print(f"Demultiplexing to {len(samples)} samples...", file=sys.stderr)

    any_output = bam_out is not None or output_pattern is not None or bam_unassigned is not None
    f_out, f_unassigned, f_samples = None, None, None
    n_keep, n_remove, n_unassigned = 0, 0, 0
    with contextlib.ExitStack() as stack:
        f_in = stack.enter_context(pysam.AlignmentFile(bam_in, "rb"))
        if any_output:
            header = f_in.header.to_dict()
            if sample_tag is not None and sample_tag == "RG":
                header = new_header_with_rg(header, samples)
            if bam_out:
                f_out = stack.enter_context(pysam.AlignmentFile(bam_out, "wb", header=header, threads=threads))
            if bam_unassigned:
                f_unassigned = stack.enter_context(
                    pysam.AlignmentFile(bam_unassigned, "wb", header=header, threads=threads)
                )
            if output_pattern:
                # file handles for all samples
                f_samples = dict()
                for sample in samples:
                    path = output_pattern.format(sample=sample)
                    f_samples[sample] = stack.enter_context(
                        pysam.AlignmentFile(path, "wb", header=header, threads=threads)
                    )

        for read in tqdm(f_in.fetch(until_eof=True), disable=not verbose):
            # pdb.set_trace()
            present = {tag_id for tag_id, tag_def in tag_defs.items() if tag_def.matches(read, delimiter=delimiter)}

            # drop reads matching any @remove rule
            if any(r.matches(present) for r in remove_rules):
                n_remove += 1
                continue

            # create unique list of samples matched by @keep rules, in order of appearance in config file
            samples_seen = set()
            matched = [
                sample for r in keep_rules
                if (sample := r.sample) not in samples_seen
                and r.matches(present)
                and not samples_seen.add(sample)
            ]
            # samples_seen.add() always returns None; use not to turn that into True, allowing unseen samples to be
            # added to samples_seen while allowing the evaluation to continue to the r.matches() call
            # see https://stackoverflow.com/a/480227

            if not matched:
                n_unassigned += 1
                if f_unassigned:
                    f_unassigned.write(read)
                continue

            n_keep += 1
            if f_out or f_samples:
                if not allow_multiple:
                    matched = matched[:1]

                if sample_tag:
                    rg = ";".join(matched)
                    read.set_tag(sample_tag, rg, value_type="Z")

                if f_out:
                    f_out.write(read)
                if f_samples:
                    for sample in matched:
                        f_samples[sample].write(read)

    return dict(keep=n_keep, remove=n_remove, unassigned=n_unassigned)

def main():
    args = parse_arguments()

    # validate arguments
    if not args.output and not args.output_pattern:
        raise ValueError("At least one of --output or --output_pattern must be provided.")
    if args.output_pattern and "{sample}" not in args.output_pattern:
        raise ValueError('--output_pattern must include the "{sample}" field.')
    if args.sample_tag and re.match(r'[A-Za-z][A-Za-z0-9]', args.sample_tag) is None:
        raise ValueError('--sample_tag must be a 2-letter SAM/BAM tag matching [A-Za-z][A-Za-z0-9].')

    counts = demultiplex_bam(
        args.config,
        args.input if args.input else sys.stdin,
        bam_out=args.output,
        bam_unassigned=args.unassigned,
        output_pattern=args.output_pattern,
        allow_multiple=args.allow_multiple,
        sample_tag=args.sample_tag,
        delimiter=args.delimiter,
        threads=args.threads,
        verbose=args.verbose,
    )
    if args.verbose:
        print("Demultiplexing complete. Read counts:", counts, file=sys.stderr)

if __name__ == "__main__":
    main()