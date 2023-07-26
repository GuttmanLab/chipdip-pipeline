import gzip
import json
import os
import re
import shutil
import shlex
import subprocess
import sys
from contextlib import AbstractContextManager


regex_fastq = re.compile(r'(\.fastq\.gz$)|(\.fq.gz$)|(\.fastq$)|(\.fq$)')
regex_bam = re.compile(r'(\.bam$)|(\.sam$)|(\.cram$)')


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


class BashRunnerWithSharedEnvironment(AbstractContextManager):
    """
    Run multiple bash scripts with persisent environment.
    Environment is stored to "env" member between runs. This can be updated
    directly to adjust the environment, or read to get variables.

    Source: https://stackoverflow.com/a/68339760
    """

    def __init__(self, env=None):
        if env is None:
            env = dict(os.environ)
        self.env = env
        self._fd_read, self._fd_write = os.pipe()

    def run(self, cmd, **opts):
        if self._fd_read is None:
            raise RuntimeError("BashRunner is already closed")
        write_env_pycode = ";".join(
            [
                "import os",
                "import json",
                f"os.write({self._fd_write}, json.dumps(dict(os.environ)).encode())",
            ]
        )
        write_env_shell_cmd = f"{sys.executable} -c '{write_env_pycode}'"
        cmd += "\n" + write_env_shell_cmd
        result = subprocess.run(
            ["bash", "-ce", cmd], pass_fds=[self._fd_write], env=self.env, **opts
        )
        self.env = json.loads(os.read(self._fd_read, 50000).decode())
        return result

    def __exit__(self, exc_type, exc_value, traceback):
        if self._fd_read:
            os.close(self._fd_read)
            os.close(self._fd_write)
            self._fd_read = None
            self._fd_write = None

    def __del__(self):
        self.__exit__(None, None, None)


def count_lines(paths, n_processes=1, cmd_unzip=None, cmd_count=None, env=None, return_raw=False):
    '''
    Args
    - paths: list-like
        Paths to read files. All read files must be of the same file type.
    - n_processes: int. default=1
        Number of parallel processes to use.
    - cmd_unzip: str. default=None
        Command to write file contents to standard output.
        If None, defaults to 'cat'.
        Example: If pigz is available, use 'unpigz -c' instead of 'gunzip -c' for a
          small speedup for gzipped files.
    - cmd_count: str. default=None
        Command that reads from standard in lines output from cmd_unzip and
        outputs an integer.
        If None, defaults to 'wc -l'.
    - env: dict. default=None
        Environment variables for shell commands.
    - return_raw: bool. default=False
        Return unprocessed output as split lines from underlying counting subprocesses.
    '''
    if cmd_unzip is None:
        cmd_unzip = 'cat'
    if cmd_count is None:
        cmd_count = 'wc -l'
    cmd_xargs = f'xargs -n 1 -P {n_processes} sh -c'
    cmd_shell = f'{cmd_unzip} "$0" | {cmd_count} | awk -v pre="$0" ' + '\'$0=pre"\t"$0\''
    cmd = shlex.split(cmd_xargs) + [cmd_shell]
    with subprocess.Popen(['ls'] + paths, stdout=subprocess.PIPE, env=env) as p1:
        with subprocess.Popen(cmd, stdin=p1.stdout, stdout=subprocess.PIPE, env=env) as p2:
            p1.stdout.close()
            out = p2.communicate()[0].decode().strip().splitlines()
    if return_raw:
        return out
    counts = {path: int(count) for path, count in map(lambda s: s.strip().rsplit('\t', 1), out)}
    return counts


def count_reads(paths, filetype=None, cmd_unzip=None, n_processes=None, agg=True, return_raw=False, **kwargs):
    '''
    Count the number of reads in FASTQ or SAM/BAM/CRAM files.

    Assumes samtools is available.

    Args
    - paths: list-like
        Paths to read files. All read files must be of the same file type.
    - filetype: str. default=None
        Type of files to count.
        - 'fastq': files ending in .fastq, .fq., .fastq.gz, or .fq.gz
        - 'bam': files ending in .sam, .bam., or .cram
        - None: automatically detect based on the first file in paths
    - cmd_unzip: str. default=None
        Command to write file contents to standard output. See count_lines().
    - n_processes: int. default=None
        Number of parallel processes to use.
        If None, defaults to 1.
    - agg: bool. default=True
        Return sum of read counts
    - return_raw: bool. default=False
        Return unprocessed output as split lines from underlying counting subprocesses.
    - **kwargs
        Additional keyword arguments to count_lines()
        - cmd_count
        - env

    Returns: depends
    - If return_raw is True: returns list(str)
    - If agg is False: returns dict(str -> int), giving the read counts for each path in paths
    - If agg is True: returns int, the sum of read counts over all files in paths
    '''
    if n_processes is None:
        n_processes = 1
    if filetype is None:
        if regex_fastq.search(paths[0]):
            assert all(map(regex_fastq.search, paths)), 'Files with differing extensions detected.'
            filetype = 'fastq'
        elif regex_bam.search(paths[0]):
            assert all(map(regex_bam.search, paths)), 'Files with differing extensions detected.'
            filetype = 'bam'
        else:
            raise ValueError((
                f'Unsupported file extension for file {paths[0]}. '
                'Must be .fastq.gz, .fq.gz, .bam, .sam., or .cram'))
    assert str(filetype).lower() in ('fastq', 'bam')
    factor = 4 if filetype == 'fastq' else 1

    if cmd_unzip is None:
        if paths[0].endswith('.gz'):
            assert all((path.endswith('.gz') for path in paths)), \
                'Files with different compressions detected.'
            if shutil.which('unpigz'):
                cmd_unzip = 'unpigz -c'
            else:
                cmd_unzip = 'gunzip -c'
        elif filetype == 'bam':
            cmd_unzip = f'samtools view -@ {int(max(n_processes / len(paths), 1))}'

    counts = count_lines(paths, cmd_unzip=cmd_unzip, return_raw=return_raw, **kwargs)
    if return_raw:
        return counts
    if factor != 1:
        counts = {path: count / factor for path, count in counts.items()}
    return sum(counts.values()) if agg else counts
