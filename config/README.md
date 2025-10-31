<!-- see https://snakemake.github.io/snakemake-workflow-catalog/docs/catalog.html and https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html for how this repository is formatted to comply with the Snakemake standardized workflow specification and the WorkflowHub standards -->

The pipeline is configured relative to the following directories:
- <a name="working-directory">working directory</a>: in decreasing order of precedence, the path specified by the `--directory` command-line parameter passed to Snakemake, the path specified by the `workdir:` directive in the Snakefile, or the directory in which Snakemake was invoked.
- <a name="workflow-directory">workflow directory</a>: directory containing the Snakefile
- <a name="input-directory">input directory</a>: directory containing the `config/` and `resources/` folders
- <a name="output-directory">output-directory</a>: directory containing the pipeline output

For a complete description of the directory structures, and for relevant workflow profile configuration settings, see the [main repository README](../README.md).

## Configuration Files

These files are located under `<input_directory>/config/`.

1. <a name="config-yaml">`config.yaml`</a>: Pipeline configuration - YAML file containing the processing settings and paths of required input files. Paths are specified relative to the [working directory](#working-directory).
   - Required? Yes. Must be provided in the [working directory](#working-directory), or specified via `--configfile <path_to_config.yaml>` when invoking Snakemake.
   - Required keys
     - `scripts_dir`: path to scripts folder in the [workflow directory](#workflow-directory)
     - `samples`: path to [`samples.json` file](#samples-json)
     - `barcode_config`: path to [barcode config file](#config-txt) (e.g., `config.txt`)
     - `bowtie2_index`: path to [Bowtie 2 genome index](#index-bt2)
     - `cutadapt_dpm`: path to [DPM sequences](#dpm-fasta)
     - `cutadapt_oligos`: path to [Antibody ID sequences](#bpm-fasta)
     - `bead_umi_length`: integer length of bead oligo UMIs
   - Optional keys: If these keys are omitted from `config.yaml` or set to `null`, then they will take on the default values indicated.
     - `output_dir` (default = `"results"`): path to create the [output directory](#output-directory) within which all intermediate and output files are placed.
     - `temp_dir` (default = `$TMPDIR` (if set) or `"/tmp"`): path to a temporary directory with sufficient free disk space, such as used by the `-T` option of [GNU sort](https://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html)
     - `barcode_format` (default = `null`): path to [barcode format file](#format-txt) (e.g., `format.txt`). If `null`, no barcode validation is performed.
     - `conda_env` (default = `"envs/chipdip.yaml"`): either a path to a conda environment YAML file ("\*.yml" or "\*.yaml") or the name of an existing conda environment. If the path to a conda environment YAML file, Snakemake will create a new conda environment within the `.snakemake` folder of the [working directory](#working-directory). [*If a relative path is used, the path is interpreted as relative to the Snakefile.*](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management)
     - `mask` (default = `null`): path to BED file of genomic regions to ignore, such as [ENCODE blacklist regions](#blacklist-bed); reads mapping to these regions are discarded. If `null`, no masking is performed.
     - `path_chrom_map` (default = `null`): path to [chromosome name map file](#chrom-map). If `null`, chromosome renaming and filtering are skipped, and the final BAM and/or bigWig files will use all chromosome names as-is from the Bowtie 2 index.
     - `deduplication_method` (default = `"RT&start&end"`): specify keys to use for chromatin reads deduplication, in addition to the cluster barcode. Alignment positions ('start' and/or 'end') and/or the DPM tag ('RT') can be combined using '&' (AND) or '|' (OR) operators, with '&' operators taking precedence.
     - `num_chunks` (default = `2`): integer between 1 and 99 giving the number of chunks to split FASTQ files from each sample into for parallel processing
     - `generate_splitbams` (default = `false`): [boolean value](https://yaml.org/type/bool.html) indicating whether to generate separate BAM files for each antibody target
     - `min_oligos` (default = `2`): integer giving the minimum count of deduplicated antibody oligo reads in a cluster for that cluster to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `proportion` and `max_size` criteria
     - `proportion` (default = `0.8`): float giving the minimum proportion of deduplicated antibody oligo reads in a cluster for that cluster to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `min_oligos` and `max_size` criteria
     - `max_size` (default = `10000`): integer giving the maximum count of deduplicated genomic DNA reads in a cluster for that cluster to be to be assigned to the corresponding antibody target; this criteria is intersected (AND) with the `proportion` and `max_size` criteria
     - `merge_samples` (default = `false`): [boolean](https://yaml.org/type/bool.html) indicating whether to merge cluster files and target-specific BAM and bigWig files across samples
     - `binsize` (default = `false`): integer specifying bigWig binsize; set to `false` to skip bigWig generation. Only relevant if generate_splitbams is `true`.
     - `bigwig_normalization` (default = `"None"`): normalization strategy for calculating coverage from reads; passed to the `--normalizeUsing` argument for the `bamCoverage` command from the deepTools suite. As of version 3.5.2, deepTools `bamCoverage` currently supports `RPKM`, `CPM`, `BPM`, `RPGC`, or `None`. Only relevant if bigWig generation is requested (i.e., `generate_splitbams` is `true` and `binsize` is not `false`).
     - `effective_genome_size` (default = `null`): integer specifying effective genome size (see [deepTools documentation](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html) for a definition). If `null`, effective genome size is computed as the number of unmasked sequences in the Bowtie 2 index, after selecting for chromosomes specified in the [chromosome name map file](#chrom-map) and excluding regions specified by the mask file. Only relevant if bigWig generation is requested using normalization strategy `RPGC` (i.e., `generate_splitbams` is `true`, `binsize` is not `false`, and `bigwig_normalization` is `RPGC`).
     - `email` (default = `null`): email to send error notifications to if errors are encountered during the pipeline. If `null`, no emails are sent.
   - Additional notes
     - `null` values can be specified explicitly (e.g., `format: null`) or implicitly (e.g., `format: `).
     - For keys `format`, `mask`, `path_chrom_map`, and `email`, an empty string `""` is treated identically to if the value is `null`.

2. <a name="samples-json">`samples.json`</a>: Samples file - JSON file with the paths of FASTQ files (read1, read2) to process.
   - Required? Yes.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `samples`
   - This can be prepared using `fastq2json.py --fastq_dir <path_to_directory_of_FASTQs>` or manually formatted as follows:

     ```json
     {
        "sample1": {
          "R1": [
            "<path_to_data>/sample1_run1_R1.fastq.gz",
            "<path_to_data>/sample1_run2_R1.fastq.gz",
          ],
          "R2": [
            "<path_to_data>/sample1_run1_R2.fastq.gz",
            "<path_to_data>/sample1_run2_R2.fastq.gz",
          ]
        },
        "sample2": {
          "R1": [
            "<path_to_data>/sample2_R1.fastq.gz"
          ],
          "R2": [
            "<path_to_data>/sample2_R2.fastq.gz"
          ]
        },
        ...
     }
     ```

   - Data assumptions:
     - FASTQ files are gzip-compressed.
     - Read names do not contain two consecutive colons (`::`). This is required because the pipeline adds `::` to the end of read names before adding barcode information; the string `::` is used as a delimiter in the pipeline to separate the original read name from the identified barcode.
   - If there are multiple FASTQ files per read orientation per sample (as shown for `sample1` in the example above), the pipeline will concatenate them and process them together as the same sample.
   - Each sample is processed independently, generating independent BAM files and statistics for quality assessment (barcode identification efficiency, cluster statistics, cluster size distributions, splitbam statistics). For ease of comparison, all samples are overlaid together in quality assessment plots.
   - The provided sample read files under the `data/` folder were simulated via a [Google Colab notebook](https://colab.research.google.com/drive/1CyjY0fJSiBl4vCz6FGFuT3IZEQR5XYlI). The genomic DNA reads correspond to ChIP-seq peaks on chromosome 19 (mm10) for transcription factors MYC (simulated as corresponding to Antibody ID `BEAD_AB1-A1`) and TCF12 (simulated as corresponding to Antibody ID `BEAD_AB2-A2`).
   - Sample names (the keys of the samples JSON file) cannot contain any periods (`.`). This is enforced to simplify wildcard pattern matching in the Snakefile and to allow the use of periods to delimit tags in a barcode string.

3. <a name="config-txt">`config.txt`</a>: Barcode config file - text file containing the sequences of split-pool tags and the expected split-pool barcode structure.
   - Required? Yes.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `barcode_config`
   - Used by: `scripts/java/BarcodeIdentification_v1.2.0.jar` (Snakefile `rule barcode_id`), `scripts/python/fastq_to_bam.py` (Snakefile `rule fastq_to_bam`), and `scripts/python/barcode_identification_efficiency.py` (Snakefile `rule barcode_identification_efficiency`).
     - This file is also parsed in the Snakefile itself to determine the length of the barcode (i.e., the number of rounds of barcoding) and if `generate_splitbams` is set to `true` in [`config.yaml`](#config-yaml), the set of antibody targets for which to generate individual de-multiplexed BAM files (and bigWig file too, if requested).
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
       # 2. Tag name: must contain "DPM", such as "DPM<xxx>"; must NOT contain "BEAD"
       #    - Can only contain alphanumeric characters, underscores, and hyphens,
       #      i.e., must match the regular expression "[a-zA-Z0-9_\-]+"
       # 3. Tag sequence (see resources/dpm96.fasta)
       # 4. Tag error tolerance: acceptable Hamming distance between
       #    expected tag sequence (column 3) and tag sequence in the read
       DPM	DPMBot6_1-A1	TGGGTGTT	0
       DPM	DPMBot6_2-A2	TGACATGT	0
       ...
       
       # Antibody ID sequences formatted as tab-delimited lines
       # - Identical format as for DPM tag sequences, except that Tag name (column 2)
       #   must start with "BEAD_".
       # - Tag sequences must match resources/bpm.fasta
       DPM	BEAD_AB1-A1	GGAACAGTT	0
       DPM	BEAD_AB2-A2	CGCCGAATT	0
       ...
       
       # Split-pool tag sequences: same 4-column tab-delimited format as the
       #   DPM and Antibody ID section above, except that
       #   Tag category (column 1) is now ODD, EVEN, or Y.
       #   Tag name must NOT contain "BEAD" or "DPM".
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
     - Because only the first 2 antibody IDs are included in the example dataset, the other antibody ID rows are commented out. This prevents generation of empty (0 byte) placeholder files for the other 94 antibody IDs.
     - Sequences
       - The design of a DPM tags allows for 9 bp of unique sequence, but only 8 bp are used in the published SPRITE tag set (in bottom tags, the 9th bp is currently a constant `'T'`). `example_config.txt` therefore only includes the unique 8 bp sequences.
       - The design of EVEN and ODD tags allows for 17 bp of unique sequence, but only 16 bp are used in the published SPRITE tag set (in bottom tags, the 17th bp is currently a constant `'T'`). `example_config.txt` further trims the 1st unique bp from the bottom tag, leaving only the middle 15 bp unique bottom sequence.
       - The design of Y (terminal) tags allows for 9-12 bp of unique sequence.
       <!-- TODO: why are the DPM sequences in the config.txt file trimmed compared to dpm96.fasta? -->

4. <a name="format-txt">`format.txt`</a>: Barcode format file - tab-delimited text file indicating which split-pool barcode tags are valid in which round of split-pool barcoding (i.e., at which positions in the barcoding string).
   - Required? No, but highly recommended.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `barcode_format`
   - Used by: `scripts/python/split_bpm_dpm.py` (Snakefile `rule split_bpm_dpm`)
   - Column 1 indicates the zero-indexed position of the barcode string where a tag can be found.
     - Term barcode tags (Y) are position `0`; the second to last round of barcoding tags are position `1`; etc. A value of `-1` in the position column indicates that the barcode tag was not used in the experiment.
   - Column 2 indicates the name of the tag. This must be the same as the name of the tag in the [barcode config file](#config-txt). If the same tag is used in multiple barcoding rounds, then it should appear multiple times in column 2, but with different values in column 1 indicating which rounds it is used in.

5. <a name="chrom-map">`chrom_map.txt`</a>: Chromosome names map - tab-delimited text file specifying which chromosomes from the Bowtie 2 index to keep and how to rename them (if at all).
   - Required? No, but necessary if using a [blacklist mask](#blacklist-bed) that uses different chromosome names than used in the Bowtie 2 index.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `path_chrom_map`
   - Used by: `scripts/python/rename_and_filter_chr.py` (Snakefile `rule rename_and_filter_chr`, `rule merge_mask`, and `rule effective_genome_size`)
   - Column 1 specifies chromosomes (following naming convention used in the index) to keep.
     - The order of chromosomes provided here is maintained in the SAM/BAM file
       header, and consequently specifies the coordinate sorting order at the
       reference sequence level.
   - Column 2 specifies new chromosome names for the corresponding chromosomes in column 1.
   - The provided `chrom-map.txt` in this repository contains examples for retaining only canonical human or mouse chromosomes (i.e., excluding alternate loci, unlocalized sequences, and unplaced sequences) and renaming them to UCSC chromosome names (i.e., `chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`) as needed. The header of provided file also includes more detailed documentation about the specific format requirements, such as allowed characters.

## Resource Files

These files are located under `<input_directory>/resources/`.

6. <a name="bpm-fasta">`bpm.fasta`</a>: FASTA file containing the sequences of Antibody IDs
   - Required? Yes.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `cutadapt_oligos`
   - Used by: `cutadapt` (Snakefile `rule cutadapt_oligo`)
   - Each sequence should be preceded by `^` to anchor the sequence during cutadapt trimming (see Snakefile `rule cutadapt_oligo`).

7. <a name="dpm-fasta">`dpm96.fasta`</a>: FASTA file containing the sequences of DPM tags
   - Required? Yes.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `cutadapt_dpm`
   - Used by: `cutadapt` (Snakefile `rule cutadapt_dpm`)
   - Each of these sequences are 10 nt long, consisting of a unique 9 nt DPM_Bottom sequences as originally designed for SPRITE (technically, only the first 8 nt are unique, and the 9th sequence is always a `T`), plus a `T` that is complementary to a 3' `A` added to a genomic DNA sequence via dA-tailing.
<!--TODO: for chromatin read 1 - we are trimming the 5' DPM, but are we trimming the 3' DPM if the read extends beyond the DNA insert sequence? -->

8. <a name="blacklist-bed">`blacklist_hg38.bed`, `blacklist_mm10.bed`</a>: blacklisted genomic regions for ChIP-seq data
   - Required? No, but highly recommended.
   - [`config.yaml`](#config-yaml) key to specify the path to this file: `mask`
   - Used by: Snakefile `rule merge_mask`, whose output is used by `rule repeat_mask` and `rule effective_genome_size`
   - For human genome release hg38, we use [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/) from ENCODE. For mouse genome release mm10, we use [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz). These BED files use UCSC chromosome names (e.g., `chr1`, `chr2`, ...). The pipeline performs chromosome name remapping (if specified) before this step.
     - Reference paper: Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. *Sci Rep*. 2019;9(1):9354. doi:[10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z)
     - Example code used to download them into the `resources/` directory:

       ```bash
       wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz |
           zcat |
           sort -V -k1,3 > "resources/blacklist_hg38.bed"

       wget -O - https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz |
           zcat |
           sort -V -k1,3 > "resources/blacklist_mm10.bed"
       ```

9. <a name="index-bt2">`index_mm10/*.bt2`, `index_hg38/*.bt2`</a>: Bowtie 2 genome index
   - Required? Yes.
   - [`config.yaml`](#config-yaml) key to specify the path to the index: `bowtie2_index`
   - Used by: Snakefile `rule bowtie2_align` and `rule effective_genome_size`
   - If you do not have an existing Bowtie 2 index, you can download [pre-built indices](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) from the Bowtie 2 developers:

     ```bash
     # for human primary assembly hg38
     mkdir -p resources/index_hg38
     wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
     unzip -j -d resources/index_hg38 GRCh38_noalt_as.zip \*.bt2

     # for mouse primary assembly mm10
     mkdir -p resources/index_mm10
     wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
     unzip -j -d resources/index_mm10 mm10.zip \*.bt2
     ```

     This will create a set of files under `resources/index_hg38` or `resources/index_mm10`. If we want to use the `mm10` genome assembly, for example, the code above will populate `resources/index_mm10` with the following files: `mm10.1.bt2`, `mm10.2.bt2`, `mm10.3.bt2`, `mm10.4.bt2`, `mm10.rev.1.bt2`, `mm10.rev.2.bt2`. The path prefix to this index (as accepted by the `bowtie2 -x <bt2-idx>` argument) is therefore `resources/index_mm10/mm10`, which is set in the configuration file, [`config.yaml`](#config-yaml).

     Note that the pre-built indices linked above use [UCSC chromosome names](https://genome.ucsc.edu/FAQ/FAQgenes.html) (`chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`). If your alignment indices use different chromosome names (e.g., Ensembl chromosome names are `1`, `2`, ..., `X`, `Y`, `MT`), update [`chrom-map.txt`](#chrom-map) such that chromosome names in BAM files are converted to UCSC chromosome names. You can check the names of the reference sequences used to build the index by using the command `bowtie2-inspect -n <bt2-idx>`.
