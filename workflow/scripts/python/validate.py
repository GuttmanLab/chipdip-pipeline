import argparse
import json
import os
import re
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import parse_chrom_map, GZIP_MAGIC_NUMBER, file_open
import yaml


REQUIRED_KEYS = (
    'scripts_dir',
    'samples',
    'barcode_config',
    'bowtie2_index',
    'cutadapt_dpm',
    'cutadapt_oligos',
    'bead_umi_length'
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Check for common mistakes in configuring ChIP-DIP pipeline"
    )
    parser.add_argument(
        "-c",
        "--config",
        metavar="FILE",
        help=("path to JSON or YAML file of pipeline config. Parameters provided via additional arguments (e.g., "
              "--chrom_map) will supersede values in the config file."),
    )
    parser.add_argument(
        "--chrom_map",
        metavar="FILE",
        help="path to chromosome name map file (e.g., chrom_map.json)",
    )
    parser.add_argument(
        "--bt2_index_summary",
        metavar="FILE",
        help="path to Bowtie 2 index summary (i.e., bowtie2-inspect --summary <index>)",
    )
    parser.add_argument(
        "--mask",
        metavar="FILE",
        help="path to mask BED file",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="suppress output"
    )
    return parser.parse_args()


def validate_samples(path_samples):
    '''
    Validate samples.json file
    '''
    with open(path_samples) as f:
        files = json.load(f)
    for sample, d in files.items():
        assert '.' not in sample, \
            f"Error in {path_samples}. Sample names are not allowed to contain periods."
        assert 'R1' in d and 'R2' in d, \
            f"Error in {path_samples}. Each sample must have read 1 (R1) and read 2 (R2) files."
        for path in d["R1"] + d["R2"]:
            with open(path, "rb") as f:
                assert f.read(2) == GZIP_MAGIC_NUMBER, \
                    f"Error in {path_samples}. FASTQ file {path} does not appear to be gzip compressed."


def validate_barcode_config(barcode_config_file):
    '''
    Validate barcode config file.
    '''
    supported_tags = ('DPM', 'EVEN', 'ODD', 'Y', 'RPM', 'LIGTAG')
    tag_layout_defined = False
    regex_tag_name = re.compile(r'[a-zA-Z0-9_\-]+')
    base_error_msg = "Error parsing line {} in the barcode config file. {}\n\t{}"
    with open(barcode_config_file, 'rt') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.upper().startswith(("#", "SPACER = ", "LAXITY = ")) or line == "":
                continue
            if line.upper().startswith("READ1 = "):
                tag_layout_defined = True
                assert all(x.upper() in supported_tags for x in line.split("= ")[1].split("|SPACER|")), \
                    base_error_msg.format(i, "Unsupported tag in READ1 tag layout.", line)
            elif line.startswith("READ2 = "):
                tag_layout_defined = True
                assert all(x.upper() in supported_tags for x in line.split("= ")[1].split("|SPACER|")), \
                    base_error_msg.format(i, "Unsupported tag in READ2 tag layout.", line)
            elif line.startswith(tuple(f"{x.upper()}\t" for x in supported_tags)):
                row = line.split("\t")
                assert len(row) == 4, \
                    base_error_msg.format(i, f"Expected 4 tab-delimited columns but only found {len(row)}.", line)
                assert regex_tag_name.match(row[1]) is not None, \
                    base_error_msg.format(i, f"Tag name did not match the pattern {regex_tag_name.pattern}.", line)
                assert re.match('[0-9]+', row[3]) is not None, \
                    base_error_msg.format(i, f"Tag error tolerance must be an integer.", line)
                if row[0].upper() == 'DPM':
                    assert ('DPM' in row[1] and (not 'BEAD' in row[1])) or row[1].startswith('BEAD_'), \
                        base_error_msg.format(
                            i,
                            (
                                "DPM tag names must contain 'DPM' and not 'BEAD', "
                                "while antibody ID tag names must start with 'BEAD_'."
                            ),
                            line
                        )
                else:
                    assert not ('DPM' in row[1] or 'BEAD' in row[1]), \
                        base_error_msg.format(i, f"Non-DPM, non-antibody ID tag names cannot contain 'BEAD' or 'DPM'.", line)
            else:
                raise ValueError(base_error_msg.format(i, "Unrecognized tag category.", line))
    assert tag_layout_defined, (
        "Barcode config file must define the tag layout for read 1 and/or read 2 orientations "
        "(i.e., lines that start with 'READ1 = ' or 'READ2 = ')."
    )


def validate_config(config):
    # check that temp_dir is writeable
    pass


def validate_format():
    pass


def validate_chrom_map(path_chrom_map, path_bt2_index_summary, chrom_map=None, chrom_sizes=None, verbose=True):
    """
    Validate the chromosome name map against the Bowtie 2 index - i.e.,
    all chromosome names to be renamed should be in the Bowtie 2 index.

    Args
    - path_chrom_map: str
        Path to chromsome name map file.
    - path_bt2_index_summary: str
        Path to Bowtie2 index summary, as generated by `bowtie2-inspect --summary <bt2_base>`
    - chrom_map: dict[str, str]. default=None
        Mapping from old chromosome name to new chromosome name.
        If not provided, generate by parsing the chromosome name map file.
    - chrom_sizes: dict[str, int]. default=None
        Mapping from chromosome name to length.
        If not provided, generate by parsing the Bowtie 2 index summary.
    - verbose: bool. default=True
        If True, print the genome size of the Bowtie 2 index and the chromosome map.

    Returns None if chromosome name map is valid.
    Raises AssertionError if chromosome name map is invalid.
    """
    if chrom_sizes is None:
        chrom_sizes = parse_bt2_index_summary(path_bt2_index_summary)
    if chrom_map is None:
        chrom_map = parse_chrom_map(path_chrom_map)
    genome_size = 0
    for chrom in chrom_map:
        assert chrom in chrom_sizes, (
            "Chromosome '{}' specified in chromosome name map '{}' not found in "
            "Bowtie 2 index summary '{}'"
        ).format(chrom, path_chrom_map, path_bt2_index_summary)
        genome_size += chrom_sizes[chrom]
    if verbose:
        print("Genome size of Bowtie 2 index:", sum(chrom_sizes.values()))
        print("Genome size of chromosome map:", genome_size)


def validate_mask(path_mask: str, chrom_map=None, chrom_sizes=None) -> None:
    """
    Check that chromosome names in the mask file match those in the chromosome name map (if provided) or
    in the Bowtie2 index (if the chromosome name map is not provided).

    At least one of chrom_map or chrom_sizes must be provided.
    If chrom_map is provided, chrom_sizes is ignored.

    Args
    - path_mask: path to BED file (may be gzip-compressed)
    - chrom_map: chromosome name map
    - chrom_sizes: chromosome sizes (obtained from Bowtie2 index summary)
      - Ignored if chrom_map is provided

    Returns: None

    Raises: AssertionError if any chromosome names in the mask file do not match those in the chromosome name map
    or in the Bowtie2 index.
    """
    if chrom_map is None and chrom_sizes is None:
        raise ValueError(
            "Either a chromosome name map or a Bowtie2 index summary must be provided to validate the mask file."
        )
    if chrom_map is not None:
        chrom_name_source = 'chromosome name map'
        valid_chroms = set(chrom_map.values())
    else:
        chrom_name_source = 'Bowtie2 index'
        valid_chroms = set(chrom_sizes.keys())

    mask_chroms_not_in_valid_chroms = set()
    at_least_one_valid_chrom = False
    with file_open(path_mask, mode='rt') as f:
        for line in f:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            chrom, *_ = line.split()
            if chrom not in valid_chroms:
                mask_chroms_not_in_valid_chroms.add(chrom)
            else:
                at_least_one_valid_chrom = True
    for chrom in mask_chroms_not_in_valid_chroms:
        print(f"Chromosome '{chrom}' in mask file '{path_mask}' not found in {chrom_name_source}.", file=sys.stderr)
    assert at_least_one_valid_chrom, \
        f'Mask BED file {path_mask} does not contain any regions from chromosomes in the {chrom_name_source}.'


def parse_bt2_index_summary(path_bt2_index_summary):
    """
    Parse the output of `bowtie2-inspect --summary <bt2_base>` to a dictionary mapping
    chromosome name to length. The chromosome name is truncated at the first whitespace,
    per SAM format specifications.
    """
    REGEX_BT2_SEQ_LINE = re.compile(r"^Sequence-\d+")
    chrom_sizes = dict()
    with open(path_bt2_index_summary, "rt") as f:
        for line in f:
            if not REGEX_BT2_SEQ_LINE.search(line):
                continue
            parts = line.strip().split("\t")
            name = parts[1]
            refname = name.split()[0]
            length = int(parts[-1])
            chrom_sizes[refname] = length
    return chrom_sizes


def main():
    args = parse_arguments()
    verbose = not args.quiet

    config = dict()
    if args.config:
        with open(args.config, "rt") as f:
            if args.config.lower().endswith(".json"):
                config = json.load(f)
            elif args.config.lower().endswith(".yaml") or args.config.lower().endswith(".yml"):
                config = yaml.safe_load(f)
            else:
                raise ValueError(("Unrecognized file extension (not .json, .yaml, or .yml) for config file: "
                                  "{}".format(args.config)))
    assert all(key in config for key in REQUIRED_KEYS), \
        'Config file must contain the following required keys: {}.'.format(', '.join(REQUIRED_KEYS))
    print(f'Config file {args.config} contains all required keys.')

    path_chrom_map = args.chrom_map if args.chrom_map else config.get("path_chrom_map", None)
    chrom_map = None
    chrom_sizes = None
    if path_chrom_map:
        chrom_map = parse_chrom_map(path_chrom_map)
    if args.bt2_index_summary:
        chrom_sizes = parse_bt2_index_summary(args.bt2_index_summary)
    if chrom_map is not None and chrom_sizes is not None:
        validate_chrom_map(path_chrom_map, args.bt2_index_summary, chrom_map=chrom_map, chrom_sizes=chrom_sizes, verbose=verbose)
        print('Validated chromosome name map file.')

    validate_barcode_config(config["barcode_config"])
    print('Validated barcode config file.')

    validate_samples(config['samples'])
    print('Validated samples file.')

    if args.mask and (chrom_map is not None or chrom_sizes is not None):
        validate_mask(config["mask"], chrom_map, chrom_sizes)
        print('Validated mask file.')
    # validate_config(config)
    # validate_format(config["format"])


if __name__ == "__main__":
    main()
