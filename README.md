Contents
- [Overview](#overview)
  - [Quick Start](#quick-start)
  - [Background](#background)
  - [Pipeline](#pipeline)
- [Input Files](#input-files)
  1. [`config.yaml`](#config-yaml)
  2. [`samples.json`](#samples-json)
  3. [`assets/bpm.fasta`](#bpm-fasta)
  4. [`assets/dpm96.fasta`](#dpm-fasta)
  5. [`config.txt`](#config-txt)
  6. [`format.txt`](#format-txt)
  7. [`assets/blacklist_hg38.bed`, `assets/blacklist_mm10.bed`](#blacklist-bed)
  8. [`assets/index_mm10/*.bt2`, `assets/index_hg38/*.bt2`](#index-bt2)
- [Output Files](#output-files)

# Overview

## Quick Start

This pipeline assumes an existing [conda](https://conda.io) installation and is written as a [Snakemake](https://snakemake.github.io/) workflow. To install Snakemake with conda, run

```
conda env create -f envs/snakemake.yaml
conda activate snakemake
```

to create and activate a conda environment named `snakemake`. Once all the [input files](#input-files) are ready, run the pipeline on a SLURM server environment with

```
./run_pipeline.sh
```

After the pipeline finishes, you can explore cluster properties in-depth using the cluster quality control Jupyter notebook (`chipdip_qc.ipynb`). The conda environment specified by `envs/chipdip_qc.yaml` includes all dependencies used by the Jupyter notebook, excluding Jupyter itself. Additionally include the dependency `jupyterlab=4.0.4` if you do not already have an existing Jupyter Notebook or Jupyter Lab setup.

Other common usage notes
- To run the pipeline for on a local computer (e.g., laptop), comment out or remove the `--cluster-config cluster.yaml` and `--cluster "sbatch ..."` arguments within `./run_pipeline.sh`, and set the number of jobs `-j <#>` to the number of local processors available.
- `run_pipeline.sh` passes any additional arguments to snakemake. For example, run `./run_pipeline.sh --dry-run` to perform a dry run, or `./run_pipeline.sh --forceall` to force (re)execution of all rules regardless of past output.
- To create a reusable ChIP-DIP conda environment, run `conda env create -f envs/chipdip.yaml` and modify the `conda_env` key in [`config.yaml`](#config-yaml) to `"chipdip"`.
- To remove intermediate pipeline files, run `./run_pipeline.sh clean`.

## Background

Terms
- **split-pool tag**: the individual sequences that are added during a single round of split-pool barcoding (DPM, EVEN, ODD, TERM)
- **split-pool barcode**: the permutation of split-pool tags that uniquely identifies a cluster
- **Antibody-ID**: a 9 nt sequence within the antibody oligo that uniquely identifies a type of antibody

<!-- TODO: figures of expected sequences -->

## Pipeline

The pipeline relies on scripts written in Java, Bash, Python. This pipeline has been validated using Java version 8.0.322 (`openjdk version "1.8.0_322"`) and Bash version 4.2.46. Versions of Python are specified in conda environments described in `envs/`, along with other third-party programs and packages that this pipeline depends on.

Workflow

0. Define samples and paths to FASTQ files (`fastq2json.py` or manually generate `samples.json`)
1. Split FASTQ files into chunks for parallel processing
2. Adaptor trimming (Trim Galore!)
3. Barcode identification 
4. Split DPM (DNA) and BPM (antibody oligo) reads into separate files
5. DPM read workflow:
   1. DPM Trimming (cutadapt)
   2. Alignment (bowtie2)
   3. Chromosome relabeling (add "chr") and filtering (removing non-canonical chromosomes)
   4. Masking (based on ENCODE blacklists)
6. BPM read workflow:
   1. BPM Trimming (cutadapt)
   2. FASTQ to BAM conversion
7. Cluster generation
8. Cluster assignment and antibody specific BAM file creation
9. Antibody-specific BigWig file creation
10. Summary plots
    1. DPM and BPM cluster size distributions
    2. Maximum representation oligo ECDFs
11. Summary statistics
    1. MultiQC (trimming, alignments)
    2. Ligation efficiency
    3. Cluster statistics
    4. Read assignment statistics

Notes / FAQ
- How are antibody oligo PCR duplicates removed? In step 6.ii above (`fastq_to_bam.py`), oligo UMIs are converted to integers and used as the reference position (column 4) in the BAM file. When clusters are generated in step 7 (during both generation of the split clusters and merged clusters), reads mapped to the same position are deduplicated.
- How are DNA sequence PCR duplicates removed? We assume that DNA molecules from a batch of cells stochastically fragmented in the crosslinking ChIP assay are unlikely to have the exact same start and end positions. Therefore, we consider identical DNA fragments (same start and end positions) in the same cluster to be PCR duplicates. As for the antibody oligo PCR duplicates, when clusters are generated in step 7 (during both generation of the split clusters and merged clusters), reads mapped to the same position are deduplicated.

See also visualizations of the pipeline generated by Snakemake (commands below assume that [Graphviz](https://graphviz.org/) is installed):
- [./dag.pdf](./dag.pdf): `snakemake --dag | dot -Tpdf > dag.pdf`
- [./filegraph.pdf](./filegraph.pdf): `snakemake --filegraph | dot -Tpdf > filegraph.pdf`
- [./rulegraph.pdf](./rulegraph.pdf): `snakemake --rulegraph | dot -Tpdf > rulegraph.pdf`

# Directory structures

We will refer to 4 directories:

1. <a name="working-directory">**Working directory**</a>: We follow the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-interpret-relative-paths) in using the term "working directory" to refer to the
   > either the directory in which you have invoked Snakemake or whatever was specified for `--directory` or the `workdir:` directive
   - This is also where Snakemake creates a `.snakemake` directory within which it installs conda environments and keeps track of metadata regarding the pipeline.

2. <a name="pipeline-directory">Pipeline directory</a>: where the software (including the `Snakefile` and scripts) resides
   - `envs/`
   - `scripts/`
   - `fastq2json.py`
   - `Snakefile`

3. <a name="input-directory">Input directory</a>: where configuration and data files reside
   - `assets/`
   - `data/`
   - `cluster.yaml`
   - [`config.yaml`](#config-yaml): paths are specified relative to the [working directory](#working-directory)
   - `config.txt` / `example_config.txt`
   - `format.txt` / `example_format.txt`
   - `samples.json` / `example_samples.json`: paths are specified relative to the [working directory](#working-directory)
   - `run_pipeline.sh`: the paths in the arguments `--snakefile <path to Snakefile>`, `--cluster-config <path to cluster.yaml>`, and `--configfile <path to config.yaml>` are relative to where you run `run_pipeline.sh`

4. <a name="output-directory">Output or workup directory</a> (`workup/`): where to place this `workup` directory can be changed in [`config.yaml`](#config-yaml)
   - `alignments/`
   - `alignments_parts/`
   - `bigwigs/`
   - `clusters/`
   - `clusters_parts/`
   - `fastqs/`
   - `logs/`
   - `qc/`
   - `splitbams/`
   - `splitfq/`
   - `trimmed/`
   - `ligation_efficiency.txt`
   - `pipeline_counts.txt`

For reproducibility, we recommend keeping the pipeline, input, and output directories together. In other words, the complete directory should look like this GitHub repository with an extra `workup` subdirectory created upon running this pipeline.

However, the pipeline directory can also be kept separate and used repeatedly on different datasets.
- The [working directory](#working-directory) should stay with the [input directory](#input-directory), so that the `.snakemake` folder containing the Snakemake pipeline metadata (that keeps track of which steps of the pipeline have completed) is paired with the configuration files.
  - Assuming that the above directory structure is followed, most of the paths in [`config.yaml`](#config-yaml) can remain relative paths to configuration files and asset files in the [input directory](#input-directory). The only paths in [`config.yaml`](#config-yaml) that need to be modified are `scripts_dir` and `output_dir`.
  - Modify the `--snakefile <path to Snakefile>` argument in `run_pipeline.sh` to point to the Snakefile in the [pipeline directory](#pipeline-directory).
  - Run `run_pipeline.sh` from the [input directory](#input-directory).
- To reuse the `chipdip` conda environment, create a discoverable `chipdip` conda environment (e.g., `conda env create -f envs/chipdip.yaml`) and set the `use_existing_conda_env` option in [`config.yaml`](#config-yaml) to `true`.

# Input Files

1. <a name="config-yaml">`config.yaml`</a>: YAML file containing the processing settings and paths of required input files. As noted [above](#input-directory), paths are specified relative to the [working directory](#working-directory).
   - `output_dir`: path to create the output directory `<output_dir>/workup` within which all intermediate and output files are placed.
   - `scripts_dir`: path to scripts folder in the [pipeline directory](#pipeline-directory)
   - `temp_dir`: path to a temporary directory, such as used by the `-T` option of [GNU sort](https://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html)
   - `bID`: path to [`config.txt` file](#config-txt)
   - `format`: path to [`format.txt` file](#format-txt)
   - `samples`: path to [`samples.json` file](#samples-json)
   - `conda_env`: either a path to a conda environment YAML file (".yml" or ".yaml") or the name of an existing conda environment. If the path to a conda environment YAML file, Snakemake will create a new conda environment within the `.snakemake` folder of the [working directory](#working-directory)
   - `cutadapt_dpm`: path to [DPM sequences](#dpm-fasta)
   - `cutadapt_oligos`: path to [antibody oligo barcode sequences](#bpm-fasta)
   - `mask`
     - `mm10`: path to mm10 genomic regions to ignore, such as [ENCODE blacklist regions](#blacklist-bed); reads mapping to these regions are discarded
     - `hg38`: path to hg38 genomic regions to ignore, such as [ENCODE blacklist regions](#blacklist-bed); reads mapping to these regions are discarded
   - `bowtie2_index`
     - `mm10`: path to [Bowtie 2 genome index](#index-bt2) for the GRCm38 (mm10) build
     - `hg38`: path to [Bowtie 2 genome index](#index-bt2) for the GRCh38 (hg38) build
   - `assembly`: currently supports either `"mm10"` or `"hg38"`
   - `num_tags`: integer giving the number of rounds of tags used, including DPM. This should equal the number of times `DPM`, `ODD`, `EVEN`, `Y` appear in the first 2 lines of the [`config.txt` file](#config-txt).
   - `num_chunks`: integer giving the number of chunks to split FASTQ files from each sample into for parallel processing
   - `generate_splitbams`: [boolean value](https://yaml.org/type/bool.html) indicating whether to generate separate BAM files for each antibody target
   - `min_oligos`: minimum count of deduplicated antibody oligo barcode reads in a cluster for that cluster to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `proportion` and `max_size` criteria
   - `proportion`: minimum proportion of deduplicated antibody oligo barcode reads in a cluster for that cluster to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `min_oligos` and `max_size` criteria
   - `max_size`: maximum count of deduplicated chromatin reads in a cluster for that cluster to be to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `proportion` and `max_size` criteria
   - `merge_samples`: [boolean value](https://yaml.org/type/bool.html) indicating whether to merge cluster files and target-specific BAM files across samples
   - `binsize`: integer specifying BigWig binsize; set to `false` to skip BigWig generation. Only relevant if generate_splitbams and merge_samples are both `true`.

2. <a name="samples-json">`samples.json`</a>: JSON file with the location of FASTQ files (read1, read2) to process.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `samples`
   - This can be prepared using `fastq2json.py --fastq_dir <path_to_directory_of_FASTQs>` or manually formatted as follows:

     ```{json}
     {
        "sample1": {
          "R1": ["<path_to_data>/sample1_R1.fastq.gz"],
          "R2": ["<path_to_data>/sample1_R2.fastq.gz"]
        },
        "sample2": {
          "R1": ["<path_to_data>/sample2_R1.fastq.gz"],
          "R2": ["<path_to_data>/sample2_R2.fastq.gz"]
        },
        ...
     }
     ```

   - The pipeline (in particular, the script `scripts/bash/split_fastq.sh`) currently only supports one read 1 (R1) and one read 2 (R2) FASTQ file per sample.
     - If there are multiple FASTQ files per read orientation per sample (for example, if the same sample was sequenced multiple times, or it was split across multiple lanes during sequencing), the FASTQ files will first need to be concatenated together, and the paths to the concatenated FASTQ files should be supplied in the JSON file.
   - Each sample is processed independently, generating independent cluster and BAM files. Statistics used for quality assessment (ligation efficiency, cluster statistics, MultiQC report, cluster size distributions, splitbam statistics) are computed independently for each sample but reported together in aggregate files to enable quick quality comparison across samples.
   - The provided sample read files under the `data/` folder were simulated via a [Google Colab notebook](https://colab.research.google.com/drive/1CyjY0fJSiBl4vCz6FGFuT3IZEQR5XYlI). The chromatin reads correspond to ChIP-seq peaks on chromosome 19 (mm10) for transcription factors MYC and TCF12.

3. <a name="bpm-fasta">`assets/bpm.fasta`</a>: FASTA file containing the sequences of antibody oligo barcodes
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `cutadapt_oligos`
   - Used by: `cutadapt` (Snakefile `rule cutadapt_oligo`)
   - Each sequence should be preceeded by `^` to anchor the sequence during cutadapt trimming (see Snakefile `rule cutadapt_oligo`).

4. <a name="dpm-fasta">`assets/dpm96.fasta`</a>: FASTA file containing the sequences of DPM tags
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `cutadapt_dpm`
   - Used by: `cutadapt` (Snakefile `rule cutadapt_dpm`)
   - Each of these sequences are 10 nt long, consisting of a unique 9 nt DPM_Bottom sequences as originally designed for SPRITE (technically, only the first 8 nt are unique, and the 9th sequence is always a `T`), plus a `T` that is complementary to a 3' `A` added to a chromatin DNA sequence via dA-tailing.
<!--TODO: for chromatin read 1 - we are trimming the 5' DPM, but are we trimming the 3' DPM if the read extends beyond the DNA insert sequence? -->

5. <a name="config-txt">`config.txt`</a>: Text file containing the sequences of split-pool tags and the split-pool barcoding setup.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `bID` (for "barcode ID")
   - Used by: `scripts/java/BarcodeIdentification_v1.2.0.jar` (Snakefile `rule barcode_id`) and `scripts/python/fastq_to_bam.py` (Snakefile `rule fastq_to_bam`)
   - Format: SPRITE configuration file (see our SPRITE [GitHub Wiki](https://github.com/GuttmanLab/sprite-pipeline/wiki/1.-Barcode-Identification#configuration-file) or [*Nature Protocols* paper](https://doi.org/10.1038/s41596-021-00633-y) for details).
     - Blank lines and lines starting with `#` are ignored.
     - An example barcoding configuration file is annotated below:

       ```
       # Barcoding layout for read 1 and read 2
       # - Y represents a terminal tag
       # - ODD, EVEN, and DPM indicate their respective tags
       # - SPACER accounts for the 7-nt sticky ends that allow ligation between tags
       READ1 = DPM
       READ2 = Y|SPACER|ODD|SPACER|EVEN|SPACER|ODD|SPACER|EVEN|SPACER|ODD
       
       # DPM tag sequences formatted as tab-delimited lines
       # 1. Tag category: DPM
       # 2. Tag name: must contain "DPM", such as "DPM<xxx>"
       # 3. Tag sequence (see assets/dpm96.fasta)
       # 4. Tag error tolerance: acceptable Hamming distance between
       #    expected tag sequence (column 3) and tag sequence in the read
       DPM	DPMBot6_1-A1	TGGGTGTT	0
       DPM	DPMBot6_2-A2	TGACATGT	0
       ...
       
       # Antibody oligo barcode sequences formatted as tab-delimited lines
       # - Identical format as for DPM tag sequences, except that the tag name (column 2)
       #   must contain "BEAD", such as "BEAD_<name of antibody>"
       DPM	BEAD_AB1-A1	GGAACAGTT	0
       DPM	BEAD_AB2-A2	CGCCGAATT	0
       ...
       
       # Split-pool tag sequences: same 4-column tab-delimited format as the 
       #   "DPM and antibody oligo barcode sequences" section above, except that 
       #   Tag category (column 1) is now ODD, EVEN, or Y
       EVEN	EvenBot_1-A1	ATACTGCGGCTGACG	2
       EVEN	EvenBot_2-A2	GTAGGTTCTGGAATC	2
       ...
       ODD	OddBot_1-A1	TTCGTGGAATCTAGC	2
       ODD	OddBot_2-A2	CCTGTGCGTTAGAGT	2
       ...
       Y	NYStgBot_1-A1	TATTATGGT	0
       Y	NYStgBot_2-A2	TAGCTACCTT	0
       ...
       ```
   - Notes regarding the entries in `example_config.txt`
     - Names: Each name ends with `#-Well` (for example, `4-A4`) where the `#` gives the row-major index of the tag in a 96-well plate, and `Well` denotes the corresponding row and column.
     - Sequences
       - Antibody ID sequences are 9 bp
       - The design of a DPM tags allows for 9 bp of unique sequence, but only 8 bp are used in the published SPRITE tag set (in bottom tags, the 9th bp is currently a constant `'T'`). `example_config.txt` therefore only includes the unique 8 bp sequences.
       - The design of EVEN and ODD tags allows for 17 bp of unique sequence, but only 16 bp are used in the published SPRITE tag set (in bottom tags, the 17th bp is currently a constant `'T'`). `example_config.txt` further trims the 1st unique bp from the bottom tag, leaving only the middle 15 bp unique bottom sequence.
       - The design of Y (terminal) tags allows for 9-12 bp of unique sequence.
       <!-- TODO: why are the DPM sequences in the config.txt file trimmed compared to dpm96.fasta? -->

6. <a name="format-txt">`format.txt`</a>: Tab-delimited text file indicating which split-pool barcode tags are valid in which round of split-pool barcoding (i.e., at which positions in the barcoding string).
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `format`
   - Used by: `scripts/python/split_dpm_bpm_fq.py` (Snakefile `rule split_bpm_dpm`)
   - Column 1 indicates the zero-indexed position of the barcode string where a tag can be found.
     - Term barcode tags (Y) are position `0`; the second to last round of barcoding tags are position `1`; etc. A value of `-1` in the position column indicates that the barcode tag was not used in the experiment.
   - Column 2 indicates the name of the tag. This must be the same as the name of the tag in [`config.txt`](#config-txt). If the same tag is used in multiple barcoding rounds, then it should appear multiple times in column 2, but with different values in column 1 indicating which rounds it is used in.

7. <a name="blacklist-bed">`assets/blacklist_hg38.bed`, `assets/blacklist_mm10.bed`</a>: blacklisted genomic regions for ChIP-seq data
   - For human genome release hg38, we use [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/) from ENCODE. For mouse genome release mm10, we use [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz).
   - Reference paper: Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep. 2019;9(1):9354. doi:10.1038/s41598-019-45839-z
   - Example code used to download them into the `assets/` directory:

     ```{bash}
     wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz |
         zcat |
         sort -V -k1,3 > "assets/blacklist_hg38.bed"

     wget -O - https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz |
         zcat |
         sort -V -k1,3 > "assets/blacklist_mm10.bed"
     ```

8. <a name="index-bt2">`assets/index_mm10/*.bt2`, `assets/index_hg38/*.bt2`</a>: Bowtie 2 genome index
   - [`config.yaml`](#config-yaml) key to specify the path to the index: `bowtie2_index: {'mm10': <mm10_index_prefix>, 'hg38': <hg38_index_prefix>}`
   - If you do not have an existing Bowtie 2 index, you can download [pre-built indices](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) from the Bowtie 2 developers:

     ```{bash}
     # for human primary assembly hg38
     mkdir -p assets/index_hg38
     wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
     unzip -j -d assets/index_hg38 GRCh38_noalt_as.zip \*.bt2

     # for mouse primary assembly mm10
     mkdir -p assets/index_mm10
     wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
     unzip -j -d assets/index_mm10 mm10.zip \*.bt2
     ```

     This will create a set of files under `assets/index_hg38` or `assets/index_mm10`. If we want to use the `mm10` genome assembly, for example, the code above will populate `assets/index_mm10` with the following files: `mm10.1.bt2`, `mm10.2.bt2`, `mm10.3.bt2`, `mm10.4.bt2`, `mm10.rev.1.bt2`, `mm10.rev.2.bt2`. The path prefix to this index (as accepted by the `bowtie2 -x <bt2-idx>` argument) is therefore `assets/index_mm10/mm10`, which is set in the configuration file, [`config.yaml`](#config-yaml).

     Note that the pre-built indices linked above use [UCSC chromosome names](https://genome.ucsc.edu/FAQ/FAQgenes.html) (`chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`). If your alignment indices use Ensembl chromosome names (`1`, `2`, ..., `X`, `Y`, `MT`), this pipeline includes a step to convert chromosome names in BAM files to UCSC chromosome names.

# Output Files

1. Barcode identification efficiency (`workup/ligation_efficiency.txt`)
   - A statistical summary of how many tags were found per read and the proportion of reads with a matching tag at each position.
   - The first type of statistic describes how many tags were identified per read. For example, consider a dataset of 10000 adapter-trimmed reads with the expected tag structure as specified in the `example_config.txt` file: 1 DPM tag on Read 1, 6 tags (Odd, Even, or Terminal) on Read 2.
     - `170 (1.7%) reads found with 1 tag.` For a small fraction reads, only 1 tag was identified; this is to be expected, whether due to ligation errors, PCR artifacts, or sequencing errors. These reads are output to `workup/fastqs/<sample_name>.part<###>.barcoded_short.fastq.gz` and are not used for analysis.
     - `7500 (75.0%) reads found with 7 tags.` This is the expected result, where all 7 tags are identified in the majority of reads. In the `split_bpm_dpm` rule, these reads are split into `workup/fastqs/<sample_name>.part<###>.barcoded_dpm.fastq.gz` or `workup/fastqs/<sample_name>.part<###>.barcoded_bpm.fastq.gz`, depending on whether a read corresponds to genomic DNA or an antibody oligo.
   - The second type of statistic describes at which positions tags were identified in the reads.
     - `9800 (98.0%) reads found with tag in position 1 (read 1, DPM).` As expected, a DPM-category tag is identified at the start of read 1.
     - `9700 (97.0%) reads found with tag in position 2 (read 2, Y).` As expected, a terminal tag is identified at the start of read 2.

2. Pipeline counts (`workup/pipeline_counts.txt`)
   - A tree-like summary of how many reads remained at each step of the pipeline, produced per aliquot and in aggregate. This can be used to quickly view the proportion of reads corresponding to antibody oligos versus chromatin reads; the proportion of properly barcoded reads; etc.
   - The 4 columns of numbers are as follows:
     1. The number of reads remaining after that step
     2. The proportion of reads remaining relative to the immediately previous step of the pipeline
     3. The proportion of reads remaining relative to the read type - antibody oligo (`bpm`) or chromatin (`dpm`)
     4. The proportion of reads remaining relative to the starting number of reads.
   - A tabular version is saved to `workup/qc/pipeline_counts.csv`

3. Cluster file (`workup/clusters/<sample>.clusters`)
   - Each line in a cluster file represents a single cluster. The first column is the cluster barcode. The remainder of the line is a tab deliminated list of reads. DNA reads are formated as `DPM[strand]_chr:start-end` and Antibody ID oligo reads are formated as `BPM[]_<AntibodyID>:<UMI>-0`.

4. Cluster statistics (`workup/clusters/cluster_statistics.txt`)
   - The number of clusters and BPM or DPM reads per library.

5. Cluster size distribtion (`workup/clusters/[BPM,DPM]_cluster_distribution.pdf`)
   - The distribution showing the proportion of clusters that belong to each size category.

6. Cluster size read distribution (`workup/clusters/[BPM,DPM]_read_distribution.pdf`)
   - The distribution showing the proportion of reads that belong to clusters of each size category. This can be more useful than the number of clusters since relatively few large clusters can contain many sequencing reads (i.e., a large fraction of the library) while many small clusters will contain few sequencing reads (i.e., a much smaller fraction of the library).
   
7. Maximum Representation Oligo ECDF (`workup/clusters/Max_representation_ecdf.pdf`)
   - A plot showing the distribution of proportion of BPM reads in each cluster that belong to the maximum represented Antibody ID in that cluster. A successful experiment should have an ECDF close to a right angle. Deviations from this indicate that beads contain mixtures of antibody ID oligos. Understanding the uniqueness of Antibody ID reads per cluster is important for choosing the thresholding parameters (`min_oligo`, `proportion`) for cluster assignment.

8. Maximum Representation Oligo Counts ECDF (`workup/clusters/Max_representation_ecdf.pdf`)
   - A plot showing the distribution of number of BPM reads in each cluster that belong to the maximum represented Antibody ID in that cluster. If clusters are nearly unique in Antibody ID composition, this plot is a surrogate for BPM size distribtuion. Understanding the number of Antibody ID reads per cluster is important for choosing the thresholding parameters (`min_oligo`, `proportion`) for cluster assignment.

9. BAM Files for Individual Antibodies (`workup/splitbams/*.bam`)
   - Thresholding criteria (`min_oligos`, `proportion`, `max_size`) for assigning individual clusters to individual antibodies are set in [`config.yaml`](#config-yaml).
   - The "none" BAM file (`<sample>.DNA.merged.labeled_none.bam`) contains DNA reads from clusters without Antibody ID reads.
   - THe "ambigious" BAM file (`<sample>.DNA.merged.labeled_ambiguous.bam`) contains DNA reads from clusters that failed the `proportion` thresholding criteria.
   - The "uncertain" BAM file (`<sample>.DNA.merged.labeled_uncertain.bam`) contains DNA reads from clusters that failed the `min_oligo` thresholding criteria.
   - The "filtered" BAM file (`<sample>.DNA.merged.labeled_filtered.bam`) contains DNA reads from clusters that failed the `max_size` thresholding criteria.

10. Read Count Summary for Individual Antibodies (`workup/splitbams/splitbam_statistics.txt`)
    - The number of read counts contained within each individual BAM file assigned to individual antibodies.

11. BigWig Files for Individual Antibodies (`workup/bigwigs/*.bw`)
    - BigWigs are generated using [`deeptools bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) with binsize set in [`config.yaml`](#config-yaml). Normalization is performed using effective genome size, which is calculated as the size of the canonical chromosomes in the Bowtie 2 index minus the size of regions in the mask for that genome.

# Credits

Adapted from the [SPRITE](https://github.com/GuttmanLab/sprite-pipeline) and [RNA-DNA SPRITE](https://github.com/GuttmanLab/sprite2.0-pipeline) pipelines by **Isabel Goronzy** ([@igoronzy](https://github.com/igoronzy)).

Updated with more extensive QC outputs and the Jupyter notebook by Benjamin Yeh  ([@bentyeh](https://github.com/bentyeh)).

Other contributors
- Andrew Perez ([@HeyDrew64](https://github.com/heydrew64))
- Mario Blanco ([@mrblanco](https://github.com/mrblanco))
- Mitchell Guttman ([@mitchguttman](https://github.com/mitchguttman))
