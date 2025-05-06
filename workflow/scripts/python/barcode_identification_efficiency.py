from collections import Counter, defaultdict
import operator
import os
import re
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import fastq_parse, file_open
import pysam

# Calculate barcode identification success rate.


def main():
    '''
    As a script, takes up to 2 arguments:
    - path to FASTQ or BAM file, where identified tags have been
        appended to read names
    - (optional) path to config.txt file
    '''
    bid = BarcodeIdentificationEfficiency()
    assert len(sys.argv) in (2, 3)
    if len(sys.argv) == 3:
        bid.get_position_names(sys.argv[2])
    bid.count_tags(sys.argv[1])
    bid.print_to_stdout()


class BarcodeIdentificationEfficiency:
    def __init__(self):
        self._aggregate_count = Counter()
        self._position_count = Counter()
        self._pattern = re.compile(r"\[([a-zA-Z0-9_\-]+)\]")
        self._total = 0
        self._position_names = defaultdict(str)

    def get_position_names(self, configfile):
        '''
        Parse read structures from config.txt file and generate
        map from tag position to position description
        '''
        tag_layout = [list(), list()]
        with open(configfile, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith(("#", "SPACER = ", "LAXITY = ")) or line == "":
                    continue
                if line.startswith("READ1 = "):
                    tag_layout[0] = line.split("= ")[1].split("|SPACER|")
                elif line.startswith("READ2 = "):
                    tag_layout[1] = line.split("= ")[1].split("|SPACER|")
                else:
                    break
        assert len(tag_layout[0] + tag_layout[1]) > 0, \
            "Could not properly parse the tag config file."
        for i, name in enumerate(tag_layout[0]):
            self._position_names[i] = f"read 1, {name}"
        for i, name in enumerate(tag_layout[1]):
            self._position_names[i + len(tag_layout[0])] = f"read 2, {name}"

    def count_tags(self, filename):
        if filename.lower().endswith(".bam"):
            self.count_tags_in_bam_file(filename)
        elif filename.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            self.count_tags_in_fastq_file(filename)
        else:
            msg = (
                f"Do not know how to open {filename}."
                "Supported extensions: .fastq, .fq, .fastq.gz, .fq.gz, .bam"
            )
            raise ValueError(msg)

    def count_tags_in_bam_file(self, bamfile):
        with pysam.AlignmentFile(bamfile, "rb") as f:
            for read in f.fetch(until_eof=True):
                self.count_tags_in_name(read.query_name)
                self._total += 1

    def count_tags_in_fastq_file(self, fastqfile):
        with file_open(fastqfile, mode="rt") as f:
            for qname, seq, thrd, qual in fastq_parse(f):
                self.count_tags_in_name(qname)
                self._total += 1

    def count_tags_in_name(self, name):
        name_split = name.split('::', 1)
        if len(name_split) == 1:
            tags = self._pattern.findall(name)
        else:
            tags = self._pattern.findall(name_split[1])
        num_found = 0
        for pos, tag in enumerate(tags):
            if tag != "NOT_FOUND":
                num_found += 1
                self._position_count[pos] += 1
        self._aggregate_count[num_found] += 1

    def print_to_stdout(self):

        counts = sorted(self._aggregate_count.items(), key=operator.itemgetter(0))

        for num_tags, count in counts:
            pct = "{0:.1f}%".format(100.0 * count / self._total)
            tag = "tag" if num_tags == 1 else "tags"
            print(f"{count} ({pct}) reads found with {num_tags} {tag}.")

        print("")
        counts = sorted(self._position_count.items(), key=operator.itemgetter(0))

        for position, count in counts:
            pct = "{0:.1f}%".format(100.0 * count / self._total)
            description = ""
            if position in self._position_names:
                description = f" ({self._position_names[position]})"
            print(
                f"{count} ({pct}) reads found with tag in position "
                f"{position + 1}{description}."
            )


if __name__ == "__main__":
    main()
