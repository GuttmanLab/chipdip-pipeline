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
