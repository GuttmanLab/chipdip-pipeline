# Overview

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