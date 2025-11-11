"""
Convenience functions related to parsing arguments and files.
"""

import argparse
import collections.abc
import contextlib
import gzip
import io
import re
import typing

# magic number for the gzip file format
# - https://en.wikipedia.org/wiki/Gzip
# - https://stackoverflow.com/a/47080739
GZIP_MAGIC_NUMBER = b"\x1f\x8b"

# valid reference sequence names, as defined in the SAM specification at
#   https://samtools.github.io/hts-specs/SAMv1.pdf
REGEX_RNAME = re.compile(r'[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*')


class AutoCloseGzipFile(gzip.GzipFile):
    """
    Subclass of gzip.GzipFile whose close() method actually closes fileobj.

    From Python documentation (https://docs.python.org/3/library/gzip.html)
    > Calling a GzipFile object’s close() method does not close fileobj,
      since you might wish to append more material after the compressed data.
    """

    def close(self):
        if self.fileobj:
            self.fileobj.close()
        super().close()


def file_open(filename: str, mode: typing.Literal["rb", "rt"] = "rb", encoding: str = "utf-8") -> typing.IO:
    """
    Detect if a file is gzip-compressed and return an appropriate file object for
    reading only (only supports "rb" and "rt" modes). Meant to be a mostly drop-in
    replacement for open()/gzip.open().
    """
    assert mode in ("rb", "rt"), 'file_open() only supports "rb" and "rt" modes'
    f = open(filename, "rb")
    first_two_bytes = f.peek(2)[:2]

    # peek is not guaranteed to return the number of bytes requested -->
    # try reading and seeking back to 0 as a fallback
    if len(first_two_bytes) != 2:
        first_two_bytes = f.read(2)
        f.seek(0)
    if first_two_bytes == GZIP_MAGIC_NUMBER:
        f = AutoCloseGzipFile(fileobj=f, mode="rb")
    if mode == "rt":
        return io.TextIOWrapper(f, encoding=encoding)
    else:
        return f


def fastq_parse(fp: typing.IO) -> collections.abc.Generator[tuple[str, str, str, str]]:
    """
    Parse FASTQ file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:
        linecount += 1
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


def fasta_parse(fp: typing.IO, header_char: str | tuple[str] = '>') -> collections.abc.Generator[tuple[str, str]]:
    """
    Parse FASTA file.

    Args
    - fp: file-like object to read from
    - header_char: character(s) indicating header lines

    Yield tuples of (header, sequence), where 'header' includes the leading header_char character.
    """
    name = None
    seq_lines = []
    for line in fp:
        line = line.strip()
        if line.startswith(header_char):
            if name is not None:
                yield name, ''.join(seq_lines)
            name = line
            seq_lines = []
        else:
            seq_lines.append(line)
    if name is not None:
        yield name, ''.join(seq_lines)


@contextlib.contextmanager
def open_path_or_file(
    file: str | typing.IO,
    opener: collections.abc.Callable[[str], typing.IO] = open,
    **kwargs
) -> collections.abc.Generator[typing.IO]:
    """
    Given either a file path or a file-like object, yield a file-like object.

    Args
    - file: if a file-like object is given, yield as-is; if a string file path is given, the opened file object is
        yielded and will be closed after use.
    - opener: function to open the file path (e.g., open, gzip.open, helpers.file_open, etc.)
    - **kwargs: keyword arguments to pass to opener()

    Yields: file-like object

    Usage:
        # open either a path or file object, transparently supporting paths to gzip-compressed files
        with open_path_or_file(file, file_open, mode="rt") as f:
            # use f here

        # pass through sys.stdin or sys.stdout as needed, without closing after use
        with open_path_or_file(sys.stdout) as f:
            f.write(...)
        # stdout is not closed after use
    """
    if isinstance(file, str):
        with opener(file, **kwargs) as f:
            yield f
    else:
        yield file


def positive_int(value: str) -> int:
    """
    Check that a string represents a positive integer.
    If so, return the integer value. Otherwise, raise an argparse.ArgumentTypeError.
    Source: https://stackoverflow.com/a/14117511
    """
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def parse_chrom_map(path: str) -> dict[str, str]:
    """
    Parse a chromosome name map file to a dict mapping old chromosome names to new names.
    """
    chrom_map = dict()
    with file_open(path, mode='rt') as f:
        for line in f:
            if line.strip() == '' or line.strip().startswith('"'):
                continue
            old_name, new_name = line.strip().split("\t")
            assert REGEX_RNAME.match(old_name) and REGEX_RNAME.match(new_name), \
                ("At least one of these chromosome names in the chromosome name map is invalid: ",
                 f"{old_name} or {new_name}")
            assert old_name not in chrom_map, \
                f"The chromosome name '{old_name}' is repeated in the chromosome name map."
            assert new_name not in chrom_map.values(), \
                f"The chromosome name '{new_name}' is repeated in the chromosome name map."
            chrom_map[old_name] = new_name
    return chrom_map
