import gzip


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
