"""
Convenience functions related to parsing arguments and files.
"""
import argparse
import gzip
import re


def file_open(filename):
    """
    Open as normal or as gzip
    """
    f = open(filename, "rb")
    if f.read(2) == b"\x1f\x8b":  # compressed always start with these two bytes
        f.seek(0)  # return to start of file
        return gzip.GzipFile(fileobj=f, mode="rb")
    else:
        f.seek(0)
        return f


def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:
        linecount += 1
        try:
            line_stripped = line.decode("UTF-8").rstrip()
        except AttributeError:
            line_stripped = line.rstrip()
        if linecount % 4 == 1:
            name = line_stripped
            assert name.startswith("@"), (
                "ERROR: The 1st line in FASTQ element does not start with '@'.\n\
                   Please check FASTQ file near line number %s"
                % (linecount)
            )
        elif linecount % 4 == 2:
            seq = line_stripped
        elif linecount % 4 == 3:
            thrd = line_stripped
            assert thrd.startswith("+"), (
                "ERROR: The 3rd line in FASTQ element does not start with '+'.\n\
                   Please check FASTQ file near line number %s"
                % (linecount)
            )
        elif linecount % 4 == 0:
            qual = line_stripped
            assert len(seq) == len(qual), (
                "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FASTQ file near line number %s"
                % (linecount)
            )
            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4


def positive_int(value):
    """
    Check that a string represents a positive integer.
    If so, return the integer value. Otherwise, raise an argparse.ArgumentTypeError.
    Source: https://stackoverflow.com/a/14117511
    """
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def parse_chrom_map(path):
    """
    Parse a chromosome name map file to a dict mapping old chromosome names to new names.
    """
    REGEX_RNAME = re.compile(r'[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*')
    chrom_map = dict()
    with open(path, 'rt') as f:
        for line in f:
            if line.strip() == '' or line.strip().startswith('"'):
                continue
            old_name, new_name = line.strip().split("\t")
            assert REGEX_RNAME.match(old_name) and REGEX_RNAME.match(new_name), \
                ("At least one of these chromosome names in the chromosome name map is invalid: ",
                 f"{old_name} or {new_name}")
            assert old_name not in chrom_map, \
                f"The chromosome name '{old_name}' is repeated in the chromosome name map."
            chrom_map[old_name] = new_name
    return chrom_map
