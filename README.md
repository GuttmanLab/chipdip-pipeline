# Overview

Terms
- **barcode**: this term is overloaded to refer to one of two possible sequences
  - A **split-pool barcode** is a concatenation of split-pool tags that uniquely identifes a cluster
  - **antibody oligo barcode**: a 9 nt sequence within the antibody oligo that uniquely identifies a type of antibody

# Required Input Files

All paths are relative to the project directory.

1. `samples.json`: JSON file with the location of FASTQ files (read1, read2) to process.
   - This can be prepared using `fastq2json.py --fastq_dir <path_to_directory_of_FASTQs>` or manually formatted as follows:
   ```{json}
   {
      "sample1": {
        "R1": ["<path_to_data>/sample1_S1_R1_001.fastq.gz"],
        "R2": ["<path_to_data>/sample1_S1_R2_001.fastq.gz"]
      },
      "sample2": {
        "R1": ["<path_to_data>/sample2_S2_R1_001.fastq.gz"],
        "R2": ["<path_to_data>/sample2_S2_R2_001.fastq.gz"]
      },
      ...
   }
   ```
   - The pipeline (in particular, the script `scripts/bash/split_fastq.sh`) currently only supports one read 1 (R1) and one read 2 (R2) FASTQ file per sample.
     - If there are multiple FASTQ files per read orientation per sample (for example, if the same sample was sequenced multiple times, or it was split across multiple lanes during sequencing), the FASTQ files will first need to be concatenated together, and the paths to the concatenated FASTQ files should be supplied in the JSON file.

2. `bpm.fasta`: FASTA file containing the sequences of antibody oligo barcodes
   - Each sequence should be preceeded by `^` to anchor the sequence during cutadapt trimming (see Snakefile `rule cutadapt_oligo`).

3. `dpm96.fasta`: FASTA file containing the sequences of DPM tags
   - Each of these sequences are 10 nt long, consisting of a unique 9 nt DPM_Bottom sequences as originally designed for SPRITE (technically, only the first 8 nt are unique, and the 9th sequence is always a `T`), plus a `T` that is complementary to a 3' `A` added to a chromatin DNA sequence via dA-tailing.
<!--TODO: for chromatin read 1 - we are trimming the 5' DPM, but are we trimming the 3' DPM if the read extends beyond the DNA insert sequence? -->

4. `config.txt`: Text file containing the sequences of split-pool tags and the split-pool barcoding setup.
   - The format of this file is identical to a SPRITE config.txt (see SPRITE nature protocols for details). 
	 - Line 1: READ1 barcoding design `READ1 = DPM`
	 - Line 2: READ2 barcoding design `READ2 = Y|SPACER|ODD|SPACER|EVEN`
	 - Lines 3+: DPM and antibody oligo barcode sequences
	   - Both DPM and antibody oligo barcodes are labeled as `DPM` in column 1. Antibody oligo barcodes should be named `BEAD_Antibody` and DPM tags should be named `DPMxxx`
     - Remaining Lines: Split-pool tags sequences `ODD`, `EVEN`, `Y`
<!--TODO: verify exactly how config.txt is parsed-->

4. `format.txt`: Tab-delimited text file indicating which split-pool barcode tags are valid in which round of split-pool barcoding (i.e., at which positions in the barcoding string).
   - Column 1 indicates the zero-indexed position of the barcode string where a tag can be found.
     - Term barcode tags (Y) are position `0`; the second to last round of barcoding tags are position `1`; etc. A value of `-1` in the position column indicates that the barcode tag was not used in the experiment.
   - Column 2 indicates the name of the tag.
   - Column 3 is the tag sequence.
   - Column 4 is the edit acceptable edit distance for this sequence.
<!--TODO: verify exactly how format.txt is parsed-->

6. `config.yaml`: YAML file containing the processing settings and locations of required input files.

# Data

## Blacklisted genomic regions for ChIP-seq data

For human genome release hg38, we use [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/) from ENCODE. For mouse genome release mm10, we use [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz).
- Reference paper: Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep. 2019;9(1):9354. doi:10.1038/s41598-019-45839-z

Example code used to download them into the `assets/` directory:

```{bash}
wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz |
    zcat |
    sort -V -k1,3 > "assets/blacklist_hg38.bed"

wget -O - https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz |
    zcat |
    sort -V -k1,3 > "assets/blacklist_mm10.bed"
```

## Indexed genomes

If you do not have an existing Bowtie 2 index, you can download [pre-built indices](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) from the Bowtie 2 developers:

```{bash}
# for human primary assembly hg38
mkdir -p assets/index_hg38
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip -d assets/index_hg38 GRCh38_noalt_as.zip \*.bt2

# for mouse primary assembly mm10
mkdir -p assets/index_mm10
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip -d assets/index_mm10 mm10.zip
```

This will create a set of files under `assets/index_hg38` or `assets/index_mm10`. If we want to use the `mm10` genome assembly, for example, the code above will populate `assets/index_mm10` with the following files: `mm10.1.bt2`, `mm10.2.bt2`, `mm10.3.bt2`, `mm10.4.bt2`, `mm10.rev.1.bt2`, `mm10.rev.2.bt2`. The path prefix to this index (as accepted by the `bowtie2 -x <bt2-idx>` argument) is therefore `assets/index_mm10/mm10`, which is set in the configuration file, `config.yaml`.

If you already have an existing Bowtie 2 index, change the corresponding value for your genome assembly under the `bowtie2_index` key in `config.yaml` to the path prefix to your index.