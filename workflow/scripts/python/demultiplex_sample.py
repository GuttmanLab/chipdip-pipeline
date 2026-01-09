"""
Demultiplex a BAM file based on sample information encoded in barcodes.

This module provides functionality to split a single BAM file into multiple BAM files (one per sample) or tag reads
with sample information based on a configurable set of rules and barcode definitions.

Configuration File Specification
--------------------------------
The configuration file is a Tab-Separated Values (TSV) file with two sections. Comments start with '#' and empty lines
are ignored.

1. Rules Section
   Specifies how to assign reads to samples or remove them.
   Directives:
     @keep <TAB> <selectors> <TAB> <sample_name>
       Assigns reads matching <selectors> to <sample_name>. <sample_name> cannot contain tab characters.
     @remove <TAB> <selectors>
       Discards reads matching <selectors>.

   Selectors:
     A comma-separated list of tag IDs or groups.
     - Tag ID: {tag_id} (e.g., {def1}) - Match specific tag definition
     - Group: {{group_name}} (e.g., {{Odd1_rowsA-C}}) - Match any tag definition in group.
     Logic: A read matches a rule if it satisfies ALL selectors in the list. A selector is satisfied if the read
            contains AT LEAST ONE of the tags specified by that selector.

2. Tag Definitions Section
   Defines the mapping between tag IDs and SAM/BAM file locations.
   Header Row (Required): id, tag, location, group (in any order; group is optional)
   Tag definition rows cannot start with the '@' character.
   Columns:
     - id: Unique identifier for the tag definition (e.g., "def1").
     - tag: The tag name that appears in the barcode string (e.g., "O1"). Multiple definitions can share the same
            tag name if they map to different locations.
     - location: Where to find the barcode in the BAM read.
                 Format: "SAM_TAG:INDEX" (e.g., "CB:1" or "RT:0").
                 INDEX is 0-based index after splitting by the delimiter.
     - group: (Optional) Group name for the tag. Use '-' for no group.

Example Config:
    @keep   {def1},{{Odd1_rowsA-C}}    sample_A
    @remove {def3}
    id      tag     location    group
    def1    DPM6    RT:0        -
    def2    O1      CB:1        Odd_rowsA-C
    def3    E1      CB:1        -


Author: Benjamin Yeh, with help from ChatGPT and Gemini 3 Pro.
"""

import argparse
import collections
import contextlib
import csv
import dataclasses
import errno
import re
import sys
import typing

import pysam
from tqdm.auto import tqdm


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Demultiplex BAM file based on encoded barcode information."
    )
    subparsers = parser.add_subparsers(dest='command', required=True, help='Subcommands')

    # --- Subcommand: demux ---
    demux_parser = subparsers.add_parser(
        'demux',
        help='Run the demultiplexing process.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    demux_parser.add_argument(
        '-i', '--input',
        help='Input BAM file. If not provided, reads from standard input.'
    )
    demux_parser.add_argument(
        '-o', '--output',
        help='Output BAM file path. At least one of --output or --output_pattern is required.'
    )
    demux_parser.add_argument(
        '--output_pattern',
        help=(
            'Format string for output BAM paths (e.g., "samples/{sample}.bam"). Must include "{sample}".'
        )
    )
    demux_parser.add_argument(
        '-u', '--unassigned',
        help='Output BAM file path for reads not matching any sample (and not removed).'
    )
    demux_parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to TSV file defining tags, groups, and rules.'
    )
    demux_parser.add_argument(
        '--allow_multiple',
        action='store_true',
        help='Allow reads to match multiple samples.'
    )
    demux_parser.add_argument(
        '--sample_tag',
        help=(
            '2-letter SAM/BAM tag to store sample info. If "RG", updates @RG header. If allow_multiple is set, values '
            'are semicolon-delimited.'
        )
    )
    demux_parser.add_argument(
        '--delimiter',
        default='.',
        help='Delimiter to split the barcode string. Empty string matches entire string.'
    )
    demux_parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for BAM compression.'
    )
    demux_parser.add_argument(
        '--progress',
        action='store_true',
        help='Show a progress bar while processing reads.'
    )

    # --- Subcommand: validate ---
    val_parser = subparsers.add_parser(
        'validate',
        help='Validate a configuration file syntax and logic.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    val_parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to TSV file defining tags, groups, and rules.'
    )

    return parser.parse_args()


# automatically add __init__, __eq__, and __hash__ methods; requires Python v3.10+
@dataclasses.dataclass(frozen=True, slots=True)
class Rule:
    """
    Represents a single keep or remove rule. A rule cannot be modified after creation.

    Attributes:
    - selectors (frozenset[frozenset[str]]): An unordered, non-redundant collection of sets of tag IDs. To match the
      rule, a read must contain at least one tag ID from EACH set in the collection.
    - sample (str | None): The sample name associated with this rule. None for @remove rules.
    """
    selectors: frozenset[frozenset[str]]
    sample: str | None

    def matches(self, present: set[str]) -> bool:
        return all(not sel.isdisjoint(present) for sel in self.selectors)


class DemuxConfig:
    """
    Parses config and maintains optimized lookup structures.

    Attributes:
    - keep_rules (list[Rule]): List of rules for assigning samples.
    - remove_rules (list[Rule]): List of rules for discarding reads.
    - lookup_map (dict): Optimized inverted index for tag lookup.
      Structure: {SAM_TAG: {INDEX: {TAG_NAME: TAG_ID}}}
      Complexity: O(K) lookup where K is number of used SAM tag/index combos.
    """

    def __init__(self, path: str):
        """
        Initializes the config parser.

        Args:
        - path: Path to the configuration TSV file.

        Raises:
        - ValueError: If the configuration file format is invalid, contains duplicates, or references undefined
          tags/groups.
        """
        self.keep_rules: list[Rule] = []
        self.remove_rules: list[Rule] = []

        # Map: SAM_TAG -> INDEX -> TAG_NAME -> TAG_ID
        self.lookup_map: dict[str, dict[int, dict[str, str]]] = \
            collections.defaultdict(lambda: collections.defaultdict(dict))

        self._parse(path)

    def _parse(self, path: str) -> None:
        """
        Parses the config file and populates rules and lookup map.

        Args:
        - path: Path to the config file.
        """
        rule_lines = []
        tag_lines = []

        try:
            with open(path) as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue
                    if s.startswith("@"):
                        rule_lines.append(s)
                    else:
                        tag_lines.append(s)
        except FileNotFoundError:
            raise ValueError(f"Config file not found: {path}")

        # 1. Parse Tag Definitions temporarily to build groups and validate existence
        all_tag_defs: dict[str, tuple[str, int, str]] = {}          # tag_id -> (sam_tag, index, tag_name)
        groups: dict[str, set[str]] = collections.defaultdict(set)  # group_name -> set of tag_ids

        reader = csv.DictReader(tag_lines, delimiter="\t")
        if not {'id', 'tag', 'location'}.issubset(reader.fieldnames or []):
            raise ValueError("Config file missing required header columns: id, tag, location")

        for row in reader:
            tag_id = row["id"]
            if tag_id in all_tag_defs:
                raise ValueError(f"Duplicate tag ID: {tag_id}")

            try:
                sam_tag, idx_str = row["location"].split(":")
                tag_name = row["tag"]
                all_tag_defs[tag_id] = (sam_tag, int(idx_str), tag_name)
            except ValueError:
                raise ValueError(f"Invalid location '{row['location']}' for tag ID '{tag_id}'. Expected 'TAG:INDEX'.")

            if (grp := row.get("group")) and grp != '-':
                groups[grp].add(tag_id)

        tag_ids_set = set(all_tag_defs.keys())

        # 2. Parse Rules with Deduplication
        seen_rules = set()

        for line in rule_lines:
            parts = line.split("\t")
            rule = None

            if line.startswith("@keep"):
                if len(parts) < 3:
                    raise ValueError(f"Invalid @keep line: {line}")
                rule = Rule(self._parse_selectors(parts[1], groups, tag_ids_set), parts[2])
                if rule not in seen_rules:
                    self.keep_rules.append(rule)
                    seen_rules.add(rule)
                else:
                    print(f"Warning: Duplicate @keep rule found: {line}", file=sys.stderr)

            elif line.startswith("@remove"):
                if len(parts) < 2:
                    raise ValueError(f"Invalid @remove line: {line}")
                rule = Rule(self._parse_selectors(parts[1], groups, tag_ids_set), sample=None)
                if rule not in seen_rules:
                    self.remove_rules.append(rule)
                    seen_rules.add(rule)
                else:
                    print(f"Warning: Duplicate @remove rule found: {line}", file=sys.stderr)

            else:
                raise ValueError(f"Unrecognized directive: {parts[0]}")

        # 3. Build optimized lookup map ONLY for used tags
        used_tag_ids = set()
        for rule in self.keep_rules + self.remove_rules:
            for selector in rule.selectors:
                used_tag_ids.update(selector)

        for tag_id in used_tag_ids:
            sam_tag, idx, tag_name = all_tag_defs[tag_id]
            self.lookup_map[sam_tag][idx][tag_name] = tag_id

    def _parse_selectors(
        self,
        expr: str,
        groups: dict[str, set[str]],
        tag_ids: typing.Collection[str]
    ) -> frozenset[frozenset[str]]:
        """
        Parses a selector string into a list of ID sets.

        Args:
        - expr: Comma-delimited string of selectors (e.g., "{id1},{{groupA}}").
        - groups: Dictionary mapping group names to sets of tag IDs.
        - tag_ids: Collection of all known/valid tag IDs.

        Returns:
        - unordered, non-redundant collection of sets of tag IDs

        Raises:
        - ValueError: If the selector syntax is invalid, references unknown groups, or references unknown tag IDs.
        """
        selectors = []
        for token in (t.strip() for t in expr.split(",") if t.strip()):
            current_selector = None
            if token.startswith("{{") and token.endswith("}}"):
                g_name = token[2:-2]
                if g_name not in groups:
                    raise ValueError(f"Undefined group referenced: {g_name}")
                current_selector = groups[g_name]
            elif token.startswith("{") and token.endswith("}"):
                tag_id = token[1:-1]
                if tag_id not in tag_ids:
                    raise ValueError(f"Undefined tag ID referenced: {tag_id}")
                current_selector = {tag_id}
            else:
                raise ValueError(f"Invalid selector format: {token}")

            # Deduplicate selectors:
            # If the exact same set of tags is already required by this rule, adding it again provides no new logic.
            if current_selector not in selectors:
                selectors.append(current_selector)
            else:
                print(f"Warning: Duplicate selector '{token}' found in rule '{expr}'. Ignoring.", file=sys.stderr)

        return frozenset(frozenset(s) for s in selectors)

    def get_present_tags(self, read: pysam.AlignedSegment, delimiter: str = ".") -> set[str]:
        """
        Identifies which defined tags are present in the read.

        Args:
        - read: The pysam read object.
        - delimiter: String delimiter to split tag values.

        Returns:
        - A set of tag IDs found in the read. Time Complexity: O(K) where K is the number of distinct SAM tag/index
          combinations used in the config.
        """
        present = set()

        # Iterate only through SAM tags known to contain used barcodes
        for sam_tag, index_map in self.lookup_map.items():
            try:
                val = read.get_tag(sam_tag)
            except KeyError:
                continue

            if delimiter:
                parts = val.split(delimiter)
                for idx, name_map in index_map.items():
                    if idx < len(parts) and (tag_id := name_map.get(parts[idx])):
                        present.add(tag_id)
            else:
                # No delimiter: index effectively 0 or whole string match
                if 0 in index_map and (tag_id := index_map[0].get(val)):
                    present.add(tag_id)

        return present


def new_header_with_rg(original_header: dict, samples: set[str]) -> dict:
    """
    Creates a new BAM header with @RG lines for the given samples.

    Args:
    - original_header: The header dictionary from the input BAM, as returned by pysam.AlignmentFile.header.to_dict().
    - samples: A set of sample names to add as Read Groups.

    Returns:
    - A new dictionary representing the BAM header with the 'RG' key updated.
    """
    new_header = original_header.copy()
    if samples:
        new_header['RG'] = [{'ID': s} for s in sorted(samples)]
    return new_header


def demultiplex_bam(
    config_path: str,
    bam_in: str | typing.IO,
    bam_out: str | typing.IO | None = None,
    bam_unassigned: str | typing.IO | None = None,
    output_pattern: str | None = None,
    allow_multiple: bool = True,
    sample_tag: str | None = None,
    delimiter: str = ".",
    threads: int = 1,
    progress: bool = False,
) -> dict[str, int]:
    """
    Demultiplexes a BAM file according to the config.

    Args:
    - config_path: Path to the TSV configuration file.
    - bam_in: Input BAM path or file-like object.
    - bam_out: Output BAM path for assigned reads (optional).
    - bam_unassigned: Output BAM path for unassigned reads (optional).
    - output_pattern: Pattern for per-sample output files (optional).
    - allow_multiple: Whether to assign reads to multiple matched samples.
    - sample_tag: SAM tag to write sample assignments into.
    - delimiter: Delimiter for splitting barcode strings.
    - threads: Number of threads for BAM compression/writing.
    - progress: If True, display a progress bar (requires tqdm).

    Returns:
    - Counts of 'keep', 'remove', and 'unassigned' reads.

    Raises:
    - ValueError: On config errors.
    - OSError: If too many output files are opened (EMFILE) or other IO errors.
    """

    cfg = DemuxConfig(config_path)
    samples = {r.sample for r in cfg.keep_rules}
    any_assigned_output = bam_out or output_pattern

    n_keep, n_remove, n_unassigned = 0, 0, 0

    with contextlib.ExitStack() as stack:
        f_in = stack.enter_context(pysam.AlignmentFile(bam_in, "rb", check_sq=False))

        f_out = None
        f_unassigned = None
        f_samples = None

        if any_assigned_output or bam_unassigned:
            header = f_in.header.to_dict()
            if sample_tag == "RG":
                header = new_header_with_rg(header, samples)

            if bam_out:
                f_out = stack.enter_context(
                    pysam.AlignmentFile(bam_out, "wb", header=header, threads=threads)
                )
            if bam_unassigned:
                f_unassigned = stack.enter_context(
                    pysam.AlignmentFile(bam_unassigned, "wb", header=header, threads=threads)
                )
            if output_pattern:
                f_samples = {}
                # Note: We do not catch OSError here anymore; we let it bubble up to main()
                for sample in samples:
                    path = output_pattern.format(sample=sample)
                    f_samples[sample] = stack.enter_context(
                        pysam.AlignmentFile(path, "wb", header=header, threads=threads)
                    )

        for read in tqdm(f_in.fetch(until_eof=True), desc="Processing reads", unit=" reads", disable=not progress):
            present = cfg.get_present_tags(read, delimiter)

            if any(r.matches(present) for r in cfg.remove_rules):
                n_remove += 1
                continue

            # Identify matching samples (preserving order, unique only)
            seen_samples = set()
            matched_samples = [
                sample
                for rule in cfg.keep_rules
                if (sample := rule.sample) not in seen_samples
                and rule.matches(present)
                and not seen_samples.add(sample)
            ]

            if not matched_samples:
                n_unassigned += 1
                if f_unassigned:
                    f_unassigned.write(read)
                continue

            n_keep += 1

            if not allow_multiple:
                matched_samples = matched_samples[:1]

            if any_assigned_output:
                if sample_tag:
                    read.set_tag(sample_tag, ";".join(matched_samples), value_type="Z")

                if f_out:
                    f_out.write(read)

                if f_samples:
                    for sample in matched_samples:
                        f_samples[sample].write(read)

    return {"keep": n_keep, "remove": n_remove, "unassigned": n_unassigned}


def main():
    """
    Main execution entry point.
    """
    args = parse_arguments()

    # --- VALIDATE COMMAND ---
    if args.command == 'validate':
        try:
            DemuxConfig(args.config)
            print(f"Validation successful: '{args.config}' is a valid config file.")
        except Exception as e:
            sys.exit(f"Validation failed: {e}")
        return

    # --- DEMUX COMMAND ---
    # Argument Validation
    if not args.output and not args.output_pattern:
        raise ValueError("At least one of --output or --output_pattern must be provided.")

    if args.output_pattern and "{sample}" not in args.output_pattern:
        raise ValueError('--output_pattern must include the "{sample}" field.')

    if args.sample_tag and not re.match(r'[A-Za-z][A-Za-z0-9]$', args.sample_tag):
        raise ValueError('--sample_tag must be a 2-letter SAM/BAM tag.')

    if args.threads < 1:
        raise ValueError("--threads must be at least 1.")

    # Use sys.stdin if input is not provided
    input_source = args.input if args.input else sys.stdin

    try:
        counts = demultiplex_bam(
            args.config,
            input_source,
            bam_out=args.output,
            bam_unassigned=args.unassigned,
            output_pattern=args.output_pattern,
            allow_multiple=args.allow_multiple,
            sample_tag=args.sample_tag,
            delimiter=args.delimiter,
            threads=args.threads,
            progress=args.progress
        )
        print(counts, file=sys.stderr)
    except ValueError as e:
        sys.exit(f"Configuration/Input Error: {e}")
    except OSError as e:
        if e.errno == errno.EMFILE:
            sys.exit(
                "\nERROR: Too many samples for system open file limit. "
                "Increase 'ulimit -n' or use --output/--sample_tag.\n"
            )
        sys.exit(f"IO Error: {e}")


if __name__ == "__main__":
    main()