"""
Aim: A Snakemake workflow to process CHIP-DIP data
"""

import json
import os
import re
import sys
import datetime
import pandas as pd

##############################################################################
# Initialize settings
##############################################################################

# Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime("%Y.%m.%d")

# Priority (lowest-to-highest) of defining configuration parameters, where each option adds to the `config` dictionary
# available in the Snakefile.
# 1. Configuration file specified by the `configfile` directive in this Snakefile: `configfile: <path_to_configfile>`
# 2. Configuration files specified on the command line: `snakemake --configfile <path_to_configfile>`
# 3. Parameters specified directly on the command line: `snakemake --config key=value`
#
# If the Snakefile includes a `configfile` directive, then a configuration file must be provided:
# - at the path specified by the `configfile` directive,
# - via the `--configfile` or `--configfiles` command line arguments,
# - or both.
# Snakemake can run the Snakefile even if the path specified by the `configfile` directive does not actually exist, as
# long as a config file is provided via the command line.
#
# See https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html.
configfile: "config.yaml"

##############################################################################
# Load required settings
##############################################################################

barcode_config = config.get("barcode_config")
if barcode_config not in (None, ""):
    print("Using barcode config:", barcode_config, file=sys.stderr)
else:
    print("Missing barcode config (barcode_config) in config.yaml", file=sys.stderr)
    sys.exit()

samples = config.get("samples")
if samples not in (None, ""):
    print("Using samples file:", samples, file=sys.stderr)
else:
    print("Missing samples file (samples) in config.yaml", file=sys.stderr)
    sys.exit()

DIR_SCRIPTS = config.get("scripts_dir")
if DIR_SCRIPTS is None:
    print("Scripts directory (scripts_dir) not specificed in config.yaml", file=sys.stderr)
    sys.exit()

def get_num_tags(path_config):
    """Parse a BarcodeID config file and return the number of tags (DPM, ODD, EVEN, Y) as an integer."""
    num_tags = 0
    with open(path_config, 'rt') as f:
        n_lines_processed = 0
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            if n_lines_processed >= 2:
                break
            line = line.strip().upper()
            if line.startswith("READ1") or line.startswith("READ2"):
                num_tags += line.count("DPM") + line.count("ODD") + line.count("EVEN") + line.count("Y")
            n_lines_processed += 1
    return num_tags

try:
    num_tags = get_num_tags(barcode_config)
    print("Using", num_tags, "tags", file=sys.stderr)
except:
    print("Could not determine number of tags from BarcodeID config file.", file=sys.stderr)
    sys.exit()

try:
    adapters = "-g file:" + config["cutadapt_dpm"]
    print("Using cutadapt sequence file", adapters, file=sys.stderr)
except:
    print("DPM adaptor sequences not specificed in config.yaml", file=sys.stderr)
    sys.exit()

try:
    oligos = "-g file:" + config["cutadapt_oligos"]
    print("Using bead oligo file", oligos, file=sys.stderr)
except:
    print("Oligo sequences not specified in config.yaml", file=sys.stderr)
    sys.exit()

bead_umi_length = config.get("bead_umi_length")
if bead_umi_length is not None:
    bead_umi_length = int(bead_umi_length)
    print("Using bead UMI length:", bead_umi_length, file=sys.stderr)
else:
    print("Bead oligo UMI length not specified in config.yaml", file=sys.stderr)
    sys.exit()

bowtie2_index = config.get("bowtie2_index")
if bowtie2_index is None:
    print("Bowtie 2 index not specified in config.yaml", file=sys.stderr)
    sys.exit()

##############################################################################
# Load optional settings
##############################################################################

email = config.get("email")
if email not in (None, ""):
    print("If any errors are encountered during the pipeline, an email will be sent to:", email, file=sys.stderr)
else:
    print("Email (email) not specified in config.yaml. Will not send email on error.", file=sys.stderr)

barcode_format = config.get("barcode_format")
if barcode_format not in (None, ""):
    print("Using barcode format file:", barcode_format, file=sys.stderr)
else:
    print(
        "(WARNING) Barcode format file not specified. The pipeline will NOT ensure barcodes are valid.",
        file=sys.stderr
    )

DIR_OUT = config.get("output_dir", "workup")
print("Using output directory:", DIR_OUT, file=sys.stderr)

temp_dir = config.get("temp_dir")
if temp_dir is not None:
    print("Using temporary directory:", temp_dir, file=sys.stderr)
else:
    temp_dir = "/central/scratch/"
    print("Defaulting to temporary directory:", temp_dir, file=sys.stderr)

num_chunks = config.get("num_chunks")
if num_chunks is not None:
    num_chunks = int(num_chunks)
    print("Splitting FASTQ files into {} chunks for parallel processing".format(num_chunks),
          file=sys.stderr)
else:
    num_chunks = 2
    print("Defaulting to 2 chunks for parallel processing", file=sys.stderr)

conda_env = config.get("conda_env")
if conda_env is None:
    conda_env = "envs/chipdip.yaml"
    print("No conda environment specified. Defaulting to envs/chipdip.yaml", file=sys.stderr)
if conda_env.strip().lower().endswith(".yaml") or conda_env.strip().lower().endswith(".yml"):
    print("Will create new conda environment from", conda_env, file=sys.stderr)
else:
    print("Using existing conda environment:", conda_env, file=sys.stderr)

mask = config.get("mask")
if mask not in (None, ""):
    print("Masking reads that align to regions in:", mask, file=sys.stderr)
else:
    mask = ""
    print("(WARNING) Mask path (mask) not specified in config.yaml, no masking will be performed.", file=sys.stderr)

path_chrom_map = config.get("path_chrom_map")
if path_chrom_map in (None, ""):
    path_chrom_map = ""
    print("Chromosome names not specified, will use all chromosomes in the Bowtie 2 index.",
          file=sys.stderr)

merge_samples = config.get("merge_samples", False)

generate_splitbams = config.get("generate_splitbams", False)
if generate_splitbams:
    min_oligos = config.get("min_oligos", 2)
    proportion = config.get("proportion", 0.8)
    max_size = config.get("max_size", 10000)
    print("Will generate BAM files for individual targets using:", file=sys.stderr)
    print("\tmin_oligos:", min_oligos, file=sys.stderr)
    print("\tproportion:", proportion, file=sys.stderr)
    print("\tmax_size:", max_size, file=sys.stderr)
else:
    print("Will not generate BAM files for individual targets.", file=sys.stderr)

binsize = config.get("binsize", False)
if binsize and not generate_splitbams:
    print("Will not generate bigWigs, because split BAMs are not being generated", file=sys.stderr)
    binsize = False

bigwig_normalization = config.get("bigwig_normalization", "None")
if binsize:
    assert bigwig_normalization in ("RPKM", "CPM", "BPM", "RPGC", "None"), (
        'bigwig_normalization config parameter must be one of '
        '"RPKM", "CPM", "BPM", "RPGC", or "None".'
    )
    print(
        "Will generate bigWig files for individual targets using normalization strategy:",
        bigwig_normalization,
        file=sys.stderr
    )

effective_genome_size = config.get("effective_genome_size", None)
compute_effective_genome_size = False
if binsize and bigwig_normalization == "RPGC":
    assert type(effective_genome_size) in (int, type(None)), \
        "effective_genome_size config parameter must be an integer or null."
    if effective_genome_size is None:
        compute_effective_genome_size = True
        print("\tWill compute effective genome size from the Bowtie 2 index.", file=sys.stderr)
    else:
        print(
            f"\tUsing user-specified effective genome size of {effective_genome_size}.",
            file=sys.stderr
        )

##############################################################################
# Location of scripts
##############################################################################

split_fastq = os.path.join(DIR_SCRIPTS, "bash", "split_fastq.sh")
barcode_id_jar = os.path.join(DIR_SCRIPTS, "java", "BarcodeIdentification_v1.2.0.jar")
barcode_identification_efficiency = os.path.join(DIR_SCRIPTS, "python", "barcode_identification_efficiency.py")
split_bpm_dpm = os.path.join(DIR_SCRIPTS, "python", "split_bpm_dpm.py")
validate = os.path.join(DIR_SCRIPTS, "python", "validate.py")
rename_and_filter_chr = os.path.join(DIR_SCRIPTS, "python", "rename_and_filter_chr.py")
make_clusters = os.path.join(DIR_SCRIPTS, "python", "make_clusters.py")
merge_clusters = os.path.join(DIR_SCRIPTS, "python", "merge_clusters.py")
fastq_to_bam = os.path.join(DIR_SCRIPTS, "python", "fastq_to_bam.py")
threshold_tag_and_split = os.path.join(DIR_SCRIPTS, "python", "threshold_tag_and_split.py")
count_unmasked_bases = os.path.join(DIR_SCRIPTS, "python", "count_unmasked_bases.py")

cluster_statistics = os.path.join(DIR_SCRIPTS, "python", "generate_cluster_statistics.py")
cluster_sizes = os.path.join(DIR_SCRIPTS, "python", "get_bead_size_distribution.py")
cluster_ecdfs = os.path.join(DIR_SCRIPTS, "python", "max_representation_ecdfs_perlib.py")

pipeline_counts = os.path.join(DIR_SCRIPTS, "python", "pipeline_counts.py")

##############################################################################
# Make output directories
##############################################################################

DIR_LOGS = os.path.join(DIR_OUT, "logs")
DIR_LOGS_CLUSTER = os.path.join(DIR_LOGS, "cluster")
os.makedirs(DIR_LOGS_CLUSTER, exist_ok=True)

##############################################################################
# Get sample files and targets
##############################################################################

with open(samples) as f:
    FILES = json.load(f)
ALL_SAMPLES = sorted(FILES.keys())

NUM_CHUNKS = [f"{i:03}" for i in range(num_chunks)]

def get_targets(barcode_config_file):
    df = pd.read_csv(
        barcode_config,
        sep="\t",
        names=["Tag", "Name", "Sequence", "Number"],
        comment="#",
        skip_blank_lines=True,
        on_bad_lines="error",
    ).dropna()
    targets = [x.replace("BEAD_", "") for x in df["Name"] if "BEAD_" in x]
    return list(set(targets))
TARGETS = get_targets(barcode_config)

print("Detected the following targets in the barcode config file:", TARGETS, file=sys.stderr)
if generate_splitbams:
    TARGETS += ['ambiguous', 'none', 'uncertain', 'filtered']
    print("  Adding 'ambiguous', 'none', 'uncertain', and 'filtered' to the list of targets.", file=sys.stderr)

# check that longest file name generated will not exceed 255 characters
# - longest file name template is {sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz
# - longest sample + target name template is {sample}.DNA.merged.labeled_{target}.bam.bai
longest_sample_name = max(map(len, ALL_SAMPLES))
assert longest_sample_name <= 214, 'Sample names must be <= 214 characters'
if generate_splitbams:
    longest_target_name = max(map(len, TARGETS))
    assert longest_sample_name + longest_target_name <= 227, \
        'Sample name + target name combined must be <= 227 characters'

##############################################################################
# Logging
##############################################################################

CONFIG = [os.path.join(DIR_LOGS, "config_" + run_date + ".json")]

BID_EFFICIENCY_ALL = [os.path.join(DIR_OUT, "barcode_identification_efficiency.txt")]

MULTI_QC = [os.path.join(DIR_OUT, "qc", "multiqc_report.html")]

LOG_VALIDATE = [os.path.join(DIR_LOGS, "validate.txt")]

##############################################################################
# Trimming
##############################################################################

SPLIT_FASTQ = expand(
    os.path.join(DIR_OUT, "split_fastq", "{sample}_{read}.part_{splitid}.fastq.gz"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM = expand(
    [os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}_val_1.fq.gz"),
     os.path.join(DIR_OUT, "trimmed", "{sample}_R2.part_{splitid}_val_2.fq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

TRIM_LOG = expand(
    os.path.join(DIR_OUT, "trimmed", "{sample}_{read}.part_{splitid}.fastq.gz_trimming_report.txt"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM_RD = expand(
    [os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz"),
     os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# Barcoding
##############################################################################

BARCODEID = expand(
    os.path.join(DIR_OUT, "fastqs", "{sample}_{read}.part_{splitid}.barcoded.fastq.gz"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

SPLIT_BPM_DPM = expand(
    [os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
     os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# DNA workup
##############################################################################

Bt2_DNA_ALIGN = expand(
    os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

CHR_DNA = expand(
    os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MASKED = expand(
    os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.masked.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_DNA = expand(
    os.path.join(DIR_OUT, "alignments", "{sample}.DNA.merged.bam"),
    sample=ALL_SAMPLES)

##############################################################################
# Bead workup
##############################################################################

FASTQ_TO_BAM = expand(
    os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.BPM.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_BEAD = expand(
    os.path.join(DIR_OUT, "alignments", "{sample}.merged.BPM.bam"),
    sample=ALL_SAMPLES)

##############################################################################
# Clustering
##############################################################################

CLUSTERS = expand(
    os.path.join(DIR_OUT, "clusters_parts", "{sample}.part_{splitid}.clusters"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

CLUSTERS_MERGED = expand(
    os.path.join(DIR_OUT, "clusters", "{sample}.clusters"),
    sample=ALL_SAMPLES)

##############################################################################
# Post Clustering
##############################################################################

CLUSTER_STATISTICS = [os.path.join(DIR_OUT, "clusters", "cluster_statistics.txt")]

CLUSTER_SIZES = [os.path.join(DIR_OUT, "clusters", "DPM_read_distribution.pdf"),
         os.path.join(DIR_OUT, "clusters", "DPM_cluster_distribution.pdf"),
         os.path.join(DIR_OUT, "clusters", "BPM_cluster_distribution.pdf"),
         os.path.join(DIR_OUT, "clusters", "BPM_read_distribution.pdf")]

ECDFS = [os.path.join(DIR_OUT, "clusters", "Max_representation_ecdf.pdf"),
         os.path.join(DIR_OUT, "clusters", "Max_representation_counts.pdf")]

LABELED_BAMS = expand(
    os.path.join(DIR_OUT, "alignments", "{sample}.DNA.merged.labeled.bam"),
    sample=ALL_SAMPLES)

SPLITBAMS = expand(
    os.path.join(DIR_OUT, "splitbams", "{sample}.DNA.merged.labeled_{target}.bam"),
    sample=ALL_SAMPLES,
    target=TARGETS)

SPLITBAMS_STATISTICS = [os.path.join(DIR_OUT, "splitbams", "splitbam_statistics.txt")]

SPLITBAMS_MERGED = expand(
    os.path.join(DIR_OUT, "splitbams", "{target}.bam"),
    target=TARGETS)

BIGWIGS = expand(
    os.path.join(DIR_OUT, "bigwigs", "{sample}.DNA.merged.labeled_{target}.bw"),
    sample=ALL_SAMPLES,
    target=TARGETS)

BIGWIGS_MERGED = expand(
    os.path.join(DIR_OUT, "bigwigs", "{target}.bw"),
    target=TARGETS)

CLUSTERS_ALL = [os.path.join(DIR_OUT, "clusters", "all.clusters")]

PIPELINE_COUNTS = [os.path.join(DIR_OUT, "pipeline_counts.txt")]

FINAL = \
    MERGE_BEAD + MERGE_DNA + CLUSTER_SIZES + ECDFS + CLUSTER_STATISTICS + MULTI_QC + \
    BID_EFFICIENCY_ALL + CONFIG + PIPELINE_COUNTS

if merge_samples:
    FINAL.extend(CLUSTERS_ALL)
if generate_splitbams:
    FINAL.extend(SPLITBAMS + SPLITBAMS_STATISTICS)
    if merge_samples:
        FINAL.extend(SPLITBAMS_MERGED)
    if binsize:
        FINAL.extend(BIGWIGS)
        if merge_samples:
            FINAL.extend(BIGWIGS_MERGED)

# ALL_OUTPUTS = \
#     SPLIT_FASTQ + TRIM + BARCODEID + SPLIT_BPM_DPM + TRIM_RD + \
#     FASTQ_TO_BAM + MERGE_BEAD + \
#     Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + \
#     CLUSTERS + CLUSTERS_MERGED + CLUSTERS_ALL + \
#     CLUSTER_STATISTICS + CLUSTER_SIZES + ECDFS + \
#     SPLITBAMS + SPLITBAMS_STATISTICS + SPLITBAMS_MERGED + \
#     BIGWIGS + BIGWIGS_MERGED + \
#     CONFIG + BID_EFFICIENCY_ALL + MULTI_QC + PIPELINE_COUNTS

##############################################################################
##############################################################################
# RULE ALL
##############################################################################
##############################################################################

rule all:
    input:
        FINAL

# Send and email if an error occurs during execution
onerror:
    if email not in (None, ""):
        shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    sample = "[^\.]+",
    target = "|".join([re.escape(x) for x in TARGETS])

# remove all output, leaving just the following in the output folder:
# - bigwigs/
# - clusters/
# - qc/
# - splitbams/
# - barcode_identification_efficiency.txt
# - pipeline_counts.txt
# - effective_genome_size.txt
rule clean:
    shell:
        '''
        for path in {DIR_OUT}/*; do
            if [[ "$path" != "{DIR_OUT}/bigwigs" ]] &&
               [[ "$path" != "{DIR_OUT}/clusters" ]] &&
               [[ "$path" != "{DIR_OUT}/qc" ]] &&
               [[ "$path" != "{DIR_OUT}/splitbams" ]] &&
               [[ "$path" != "{DIR_OUT}/barcode_identification_efficiency.txt" ]] &&
               [[ "$path" != "{DIR_OUT}/effective_genome_size.txt" ]] &&
               [[ "$path" != "{DIR_OUT}/pipeline_counts.txt" ]]; then
                echo "Removing $path" && rm -rf "$path"
            fi
        done
        '''

# Output all snakemake configuration parameters into logs folder with run date
rule log_config:
    output:
        CONFIG
    run:
        with open(output[0], "wt") as f:
            json.dump(config, f, indent=4, sort_keys=True)

# Check that configuration files and assets are set up correctly
rule validate:
    input:
        config = CONFIG
    log:
        log = LOG_VALIDATE,
        bt2_sum = os.path.join(DIR_LOGS, "bowtie2_index_summary.txt"),
    conda:
        conda_env
    shell:
        '''
        {{
            bowtie2-inspect --summary "{bowtie2_index}" > "{log.bt2_sum}"
            python "{validate}" -c "{input.config}" --bt2_index_summary "{log.bt2_sum}"
        }} &> "{log.log}"
        '''

##############################################################################
# Trimming and barcode identification
##############################################################################

# Split fastq files into chunks to processes in parallel
rule split_fastq:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        temp(expand(
            [os.path.join(DIR_OUT, "split_fastq", "{{sample}}_R1.part_{splitid}.fastq"),
             os.path.join(DIR_OUT, "split_fastq", "{{sample}}_R2.part_{splitid}.fastq")],
             splitid=NUM_CHUNKS))
    log:
        os.path.join(DIR_LOGS, "{sample}.split_fastq.log")
    params:
        dir = os.path.join(DIR_OUT, "split_fastq"),
        prefix_r1 = "{sample}_R1.part_0",
        prefix_r2 = "{sample}_R2.part_0"
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            bash "{split_fastq}" {num_chunks} "{params.dir}" "{params.prefix_r1}" {threads} {input.r1:q}
            bash "{split_fastq}" {num_chunks} "{params.dir}" "{params.prefix_r2}" {threads} {input.r2:q}
        }} &> "{log}"
        '''

# Compress the split fastq files
rule compress_fastq:
    input:
        r1 = os.path.join(DIR_OUT, "split_fastq", "{sample}_R1.part_{splitid}.fastq"),
        r2 = os.path.join(DIR_OUT, "split_fastq", "{sample}_R2.part_{splitid}.fastq")
    output:
        r1 = os.path.join(DIR_OUT, "split_fastq", "{sample}_R1.part_{splitid}.fastq.gz"),
        r2 = os.path.join(DIR_OUT, "split_fastq", "{sample}_R2.part_{splitid}.fastq.gz")
    conda:
        conda_env
    threads:
        8
    shell:
        '''
        pigz -p {threads} "{input.r1}" "{input.r2}"
        '''

# Trim adaptors
rule adaptor_trimming:
    input:
        [os.path.join(DIR_OUT, "split_fastq", "{sample}_R1.part_{splitid}.fastq.gz"),
         os.path.join(DIR_OUT, "split_fastq", "{sample}_R2.part_{splitid}.fastq.gz")]
    output:
         os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}_val_1.fq.gz"),
         os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.fastq.gz_trimming_report.txt"),
         os.path.join(DIR_OUT, "trimmed", "{sample}_R2.part_{splitid}_val_2.fq.gz"),
         os.path.join(DIR_OUT, "trimmed", "{sample}_R2.part_{splitid}.fastq.gz_trimming_report.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.adaptor_trimming.log")
    params:
        dir = os.path.join(DIR_OUT, "trimmed")
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        if [[ {threads} -gt 8 ]]; then
            cores=2
        else
            cores=1
        fi

        trim_galore \\
          --paired \\
          --gzip \\
          --cores $cores \\
          --quality 20 \\
          --fastqc \\
          -o "{params.dir}" \\
          {input:q} &> "{log}"
        '''

# Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}_val_1.fq.gz"),
        r2 = os.path.join(DIR_OUT, "trimmed", "{sample}_R2.part_{splitid}_val_2.fq.gz")
    output:
        r1_barcoded = os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded.fastq.gz"),
        r2_barcoded = os.path.join(DIR_OUT, "fastqs", "{sample}_R2.part_{splitid}.barcoded.fastq.gz")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.barcode_id.log")
    conda:
        conda_env
    shell:
        '''
        java -jar "{barcode_id_jar}" \\
          --input1 "{input.r1}" --input2 "{input.r2}" \\
          --output1 "{output.r1_barcoded}" --output2 "{output.r2_barcoded}" \\
          --config "{barcode_config}" &> "{log}"
        '''

# Calculate barcode identification efficiency
rule barcode_identification_efficiency:
    input:
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        temp(os.path.join(DIR_OUT, "{sample}.part_{splitid}.bid_efficiency.txt"))
    conda:
        conda_env
    shell:
        '''
        python "{barcode_identification_efficiency}" "{input}" "{barcode_config}" > "{output}"
        '''

rule cat_barcode_identification_efficiency:
    input:
        expand(
            os.path.join(DIR_OUT, "{sample}.part_{splitid}.bid_efficiency.txt"),
            sample=ALL_SAMPLES,
            splitid=NUM_CHUNKS)
    output:
        BID_EFFICIENCY_ALL
    shell:
        '''
        tail -n +1 {input:q} > "{output}"
        '''

# Split barcoded reads into BPM and DPM, remove incomplete barcodes
rule split_bpm_dpm:
    input:
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"),
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_other.fastq.gz"),
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_short.fastq.gz")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.split_bpm_dpm.log")
    params:
        format = f"--format '{barcode_format}'" if barcode_format not in (None, "") else ""
    conda:
       conda_env
    shell:
        '''
        python "{split_bpm_dpm}" --r1 "{input}" {params.format} &> "{log}"
        '''

##############################################################################
# Cutadapt
##############################################################################

# Trim DPM from read1 of DPM reads, remove DPM dimer reads
rule cutadapt_dpm:
    input:
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz")
    output:
        fastq = os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz"),
        qc = os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.qc.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.DPM.cutadapt.log")
    params:
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA " + adapters,
        others = "--minimum-length 20"
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        {{
            cutadapt \\
              {params.adapters_r1} \\
              {params.others} \\
              -o "{output.fastq}" \\
              -j {threads} \\
              "{input}" > "{output.qc}"

            fastqc "{output.fastq}"
        }} &> "{log}"
        '''

# Trim 9mer oligo sequence from read1 of BPM reads
rule cutadapt_oligo:
    input:
        os.path.join(DIR_OUT, "fastqs", "{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz")
    output:
        fastq = os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"),
        qc = os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.BPM.cutadapt.log")
    params:
        adapters_r1 = oligos
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        cutadapt \\
          {params.adapters_r1} \\
          -o "{output.fastq}" \\
          -j {threads} \\
          "{input}" > "{output.qc}" 2> "{log}"
        '''

##############################################################################
# DNA alignment
##############################################################################

# Align DPM reads
rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    '''
    input:
        os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz")
    output:
        os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.bowtie2.log")
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        {{
            bowtie2 \\
              -p {threads} \\
              -t \\
              --phred33 \\
              -x "{bowtie2_index}" \\
              -U "{input}" |
            samtools view -bq 20 -F 4 -F 256 - |
            samtools sort -@ {threads} -o "{output}"
        }} &> "{log}"
        '''

# Rename chromosome names and filter for chromosomes of interest
rule rename_and_filter_chr:
    input:
        os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam")
    output:
        os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.rename_and_filter_chr.log")
    params:
        chrom_map = f"--chrom_map '{path_chrom_map}'" if path_chrom_map != "" else ""
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        python "{rename_and_filter_chr}" {params.chrom_map} -t {threads} --try-symlink \\
            -o "{output}" "{input}" &> "{log}"
        '''

# Merge mask
# - Merging overlapping regions should increase the speed of running bedtools intersect in the repeat_mask rule.
# - The merged mask is also used in the generate_bigwigs rule.
# - The merged mask is sorted as follows:
#   - If a chromosome name map is provided, then it is sorted by the new chromosome names (i.e., according to names in
#     the second column of the chromosome name map file). Entries with chromosome names not in the chromosome name map
#     are discarded.
#   - If a chromosome name map is not provided, then it is sorted by the order of chromosomes in the Bowtie 2 index.
rule merge_mask:
    output:
        bed = temp(os.path.join(DIR_OUT, "mask_merge.bed")),
        genome = temp(os.path.join(DIR_OUT, "mask_merge.genome"))
    log:
        os.path.join(DIR_LOGS, "merge_mask.log")
    conda:
        conda_env
    shell:
        '''
        {{
            if [ -n "{mask}" ]; then
                if [ -n "{path_chrom_map}" ]; then
                    # chromosome name map is provided --> sort chromosomes by the new chromosome names
                    sort -k1,1 -k2,2n "{mask}" |
                    bedtools merge |
                    python "{rename_and_filter_chr}" --bed --chrom_map "{path_chrom_map}" - > "{output.bed}"

                    # create genome file for bedtools intersect
                    join -t $'\\t' \\
                        <(grep -E -e '^[^"]\\s*\\S+\\s*' "{path_chrom_map}" ) \\
                        <(bowtie2-inspect -s "{bowtie2_index}" |
                        grep -E -e '^Sequence-[0-9]+' |
                        cut -f 2,3) |
                    cut -f 2,3 > "{output.genome}"
                else
                    # chromosome name map is not provided --> sort chromosomes by their order in the Bowtie 2 index
                    sort -k1,1 -k2,2n "{mask}" |
                    bedtools merge |
                    python "{rename_and_filter_chr}" \\
                        --bed \\
                        --chrom_map <(bowtie2-inspect -n "{bowtie2_index}" | sed -E 's/(\\S+)/\\1\\t\\1/') \\
                        - > "{output.bed}"

                    # create genome file for bedtools intersect
                    bowtie2-inspect -s "{bowtie2_index}" |
                    grep -E -e '^Sequence-[0-9]+' |
                    cut -f 2,3 > "{output.genome}"
                fi
            else
                touch "{output.bed}" "{output.genome}"
            fi
        }} &> "{log}"
        '''

# Repeat mask aligned DNA reads
rule repeat_mask:
    input:
        bam = os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.bam"),
        mask = os.path.join(DIR_OUT, "mask_merge.bed"),
        genome = os.path.join(DIR_OUT, "mask_merge.genome")
    output:
        os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.masked.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.repeat_mask.log")
    conda:
        conda_env
    shell:
        '''
        {{
            if [ -n "{mask}" ]; then
                # -v: only report entries in A that have no overlap in B
                bedtools intersect \\
                    -v \\
                    -a "{input.bam}" \\
                    -b "{input.mask}" \\
                    -sorted \\
                    -g "{input.genome}" > "{output}"
            else
                echo "No mask file specified, skipping masking."
                ln -s "{input.bam}" "{output}"
            fi
        }} &> "{log}"
        '''

# Combine all mapped DNA reads into a single bam file per sample
rule merge_dna:
    input:
        expand(
            os.path.join(DIR_OUT, "alignments_parts", "{{sample}}.part_{splitid}.DNA.chr.masked.bam"),
            splitid=NUM_CHUNKS)
    output:
        os.path.join(DIR_OUT, "alignments", "{sample}.DNA.merged.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_DNA.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        samtools merge -@ {threads} "{output}" {input:q} &> "{log}"
        '''

##############################################################################
# Workup Bead Oligo
##############################################################################

# Convert the BPM FASTQ reads into a BAM file, keeping the UMI
rule fastq_to_bam:
    input:
        os.path.join(DIR_OUT, "trimmed", "{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")
    output:
        sorted = os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.BPM.bam"),
        bam = temp(os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.BPM.unsorted.bam"))
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.fastq_to_bam.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            python "{fastq_to_bam}" "{input}" "{output.bam}" "{barcode_config}" "{bead_umi_length}"
            samtools sort -@ {threads} -o "{output.sorted}" "{output.bam}"
        }} &> "{log}"
        '''

# Combine all oligo reads into a single file per sample
rule merge_beads:
    input:
        expand(
            os.path.join(DIR_OUT, "alignments_parts", "{{sample}}.part_{splitid}.BPM.bam"),
            splitid=NUM_CHUNKS)
    output:
        os.path.join(DIR_OUT, "alignments", "{sample}.merged.BPM.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_beads.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        samtools merge -@ {threads} "{output}" {input:q} &> "{log}"
        '''

##############################################################################
# Make clusters
##############################################################################

# Make clusters from aligned DNA reads and oligo reads
rule make_clusters:
    input:
        dpm = os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.DNA.chr.masked.bam"),
        bpm = os.path.join(DIR_OUT, "alignments_parts", "{sample}.part_{splitid}.BPM.bam")
    output:
        unsorted = temp(os.path.join(DIR_OUT, "clusters_parts", "{sample}.part_{splitid}.unsorted.clusters")),
        sorted = os.path.join(DIR_OUT, "clusters_parts", "{sample}.part_{splitid}.clusters")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.make_clusters.log")
    threads:
        4
    conda:
        conda_env
    shell:
        '''
        {{
            python "{make_clusters}" \\
              -i "{input.bpm}" "{input.dpm}" \\
              -o "{output.unsorted}" \\
              -n {num_tags}

            sort -k 1 -T "{temp_dir}" --parallel={threads} "{output.unsorted}" > "{output.sorted}"
        }} &> "{log}"
        '''

# Merge clusters from parallel processing into a single cluster file per sample
rule merge_clusters:
    input:
        expand(
            os.path.join(DIR_OUT, "clusters_parts", "{{sample}}.part_{splitid}.clusters"),
            splitid=NUM_CHUNKS)
    output:
        mega = temp(os.path.join(DIR_OUT, "clusters", "{sample}.duplicated.clusters")),
        final = os.path.join(DIR_OUT, "clusters", "{sample}.clusters")
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_clusters.log")
    conda:
       conda_env
    threads:
        10
    shell:
        '''
        {{
            sort -k 1 -T "{temp_dir}" --parallel={threads} -m {input:q} > "{output.mega}"
            python "{merge_clusters}" -i "{output.mega}" -o "{output.final}"
        }} &> "{log}"
        '''

rule clusters_all:
    input:
        CLUSTERS_MERGED
    output:
        CLUSTERS_ALL
    log:
        os.path.join(DIR_LOGS, "clusters_all.log")
    threads:
        10
    shell:
        '''
        sort -k 1 -T "{temp_dir}" --parallel={threads} -m {input:q} > "{output}"
        '''

##############################################################################
# Profile clusters
##############################################################################

# Generate simple statistics for clusters
rule generate_cluster_statistics:
    input:
        (CLUSTERS_ALL + CLUSTERS_MERGED) if merge_samples else CLUSTERS_MERGED
    output:
        CLUSTER_STATISTICS
    log:
        os.path.join(DIR_LOGS, "cluster_statistics.log")
    params:
        dir = os.path.join(DIR_OUT, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{cluster_statistics}" --directory "{params.dir}" --pattern .clusters \\
            > "{output}" 2> "{log}"
        '''

# Generate ecdfs of oligo distribution
rule generate_cluster_ecdfs:
    input:
        (CLUSTERS_ALL + CLUSTERS_MERGED) if merge_samples else CLUSTERS_MERGED
    output:
        ecdf = os.path.join(DIR_OUT, "clusters", "Max_representation_ecdf.pdf"),
        counts = os.path.join(DIR_OUT, "clusters", "Max_representation_counts.pdf")
    log:
        os.path.join(DIR_LOGS, "cluster_ecdfs.log")
    params:
        dir = os.path.join(DIR_OUT, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{cluster_ecdfs}" --directory "{params.dir}" --pattern .clusters \\
            --xlim 30 &> "{log}"
        '''

# Profile size distribution of clusters
rule get_size_distribution:
    input:
        (CLUSTERS_ALL + CLUSTERS_MERGED) if merge_samples else CLUSTERS_MERGED
    output:
        dpm = os.path.join(DIR_OUT, "clusters", "DPM_read_distribution.pdf"),
        dpm2 = os.path.join(DIR_OUT, "clusters", "DPM_cluster_distribution.pdf"),
        bpm = os.path.join(DIR_OUT, "clusters", "BPM_read_distribution.pdf"),
        bpm2 = os.path.join(DIR_OUT, "clusters", "BPM_cluster_distribution.pdf")
    log:
        os.path.join(DIR_LOGS, "size_distribution.log")
    params:
        dir = os.path.join(DIR_OUT, "clusters")
    conda:
        conda_env
    shell:
        '''
        {{
            python "{cluster_sizes}" --directory "{params.dir}" --pattern .clusters \\
              --readtype BPM
            python "{cluster_sizes}" --directory "{params.dir}" --pattern .clusters \\
              --readtype DPM
        }} &> "{log}"
        '''

##############################################################################
# Logging and MultiQC
##############################################################################

# Aggregate metrics using multiqc
rule multiqc:
    input:
        CLUSTERS_MERGED
    output:
        os.path.join(DIR_OUT, "qc", "multiqc_report.html")
    log:
        os.path.join(DIR_LOGS, "multiqc.log")
    params:
        dir_qc = os.path.join(DIR_OUT, "qc")
    conda:
        conda_env
    shell:
        '''
        multiqc -f -o "{params.dir_qc}" "{DIR_OUT}" &> "{log}"
        '''

rule pipeline_counts:
    input:
        SPLIT_FASTQ + TRIM + BARCODEID + SPLIT_BPM_DPM + TRIM_RD + \
        FASTQ_TO_BAM + MERGE_BEAD + \
        Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + CLUSTERS + \
        CLUSTERS_MERGED + (SPLITBAMS if generate_splitbams else [])
    output:
        csv = os.path.join(DIR_OUT, "qc", "pipeline_counts.csv"),
        pretty = os.path.join(DIR_OUT, "pipeline_counts.txt")
    log:
        os.path.join(DIR_LOGS, "pipeline_counts.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        {{
            python "{pipeline_counts}" \\
              --samples "{samples}" \\
              -w "{DIR_OUT}" \\
              -o "{output.csv}" \\
              -n {threads} |
            column -t -s $'\\t' > "{output.pretty}"
        }} &> "{log}"
        '''

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule threshold_tag_and_split:
    input:
        bam = os.path.join(DIR_OUT, "alignments", "{sample}.DNA.merged.bam"),
        clusters = os.path.join(DIR_OUT, "clusters", "{sample}.clusters")
    output:
        expand(
            os.path.join(DIR_OUT, "splitbams", "{{sample}}.DNA.merged.labeled_{target}.bam"),
            target=TARGETS
        )
    log:
        os.path.join(DIR_LOGS, "{sample}.splitbams.log")
    params:
        dir_splitbams = os.path.join(DIR_OUT, "splitbams"),
        labeled_bam = os.path.join(DIR_OUT, "alignments", "{sample}.DNA.merged.labeled.bam")
    conda:
        conda_env
    shell:
        '''
        {{
            python "{threshold_tag_and_split}" \\
                "{input.bam}" \\
                "{params.labeled_bam}" \\
                "{input.clusters}" \\
                "{params.dir_splitbams}" \\
                --min_oligos {min_oligos} \\
                --proportion {proportion} \\
                --max_size {max_size} \\
                --num_tags {num_tags}

            # touch BAM files for targets that did not get any reads passing our filters
            paths=({output:q})
            for path in "${{paths[@]}}"; do
                if [ ! -f "$path" ]; then
                    target="$(echo "$path" | sed -E -e 's/.*\\.DNA\\.merged.\\labeled_(.*)\\.bam/\\1/')"
                    echo "No split BAM file generated for target ${{target}} at ${{path}}. Touching ${{path}}."
                    touch "$path"
                fi
            done
        }} &> "{log}"
        '''

# for each target, merge splitbams across samples
rule splitbams_merged:
    input:
        expand(
            os.path.join(DIR_OUT, "splitbams", "{sample}.DNA.merged.labeled_{{target}}.bam"),
            sample=ALL_SAMPLES)
    output:
        os.path.join(DIR_OUT, "splitbams", "{target}.bam")
    log:
        os.path.join(DIR_LOGS, "splitbam_{target}.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            # check which input BAM files are non-empty
            all_paths=({input:q})
            non_empty_paths=()
            for path in "${{all_paths[@]}}"; do
                if [ -s "$path" ]; then
                    non_empty_paths+=("$path")
                fi
            done

            if [ ${{#non_empty_paths[@]}} -eq 0 ]; then
                echo "All BAM files for the target are empty. Creating empty merged BAM file."
                touch "{output}"
            else
                samtools merge -f -@ {threads} "{output}" "${{non_empty_paths[@]}}"
            fi
        }} &> "{log}"
        '''

rule index_splitbams:
    input:
        os.path.join(DIR_OUT, "splitbams", "{file}.bam")
    output:
        os.path.join(DIR_OUT, "splitbams", "{file}.bam.bai")
    log:
        os.path.join(DIR_LOGS, "splitbam-index_{file}.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        if [ -s "{input}" ]; then
            # file is non-empty; assume to be valid SAM/BAM file
            samtools index -@ {threads} "{input}"
        else
            # BAM file is empty
            touch "{output}"
        fi
        '''

# Generate summary statistics of individiual bam files
rule generate_splitbam_statistics:
    input:
        SPLITBAMS + (SPLITBAMS_MERGED if merge_samples else [])
    output:
        SPLITBAMS_STATISTICS
    log:
        os.path.join(DIR_LOGS, "splitbam_statistics.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            paths=({input:q})
            for path in "${{paths[@]}}"; do
                if [ -s "$path" ]; then
                    # file is non-empty; assume to be valid SAM/BAM file
                    count=$(samtools view -@ {threads} -c "$path")
                else
                    # file is empty
                    count=0
                fi
                echo -e "${{path}}\\t${{count}}" >> "{output}"
            done
        }} &> "{log}"
        '''

# Calculate effective genome size as defined by deepTools: the number of unmasked bases in the genome.
# See https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
rule effective_genome_size:
    input:
        mask = os.path.join(DIR_OUT, "mask_merge.bed")
    output:
        os.path.join(DIR_OUT, "effective_genome_size.txt")
    log:
        os.path.join(DIR_LOGS, "effective_genome_size.log")
    params:
        chrom_map = f"--chrom_map '{path_chrom_map}'" if path_chrom_map != "" else ""
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            if [[ "{compute_effective_genome_size}" = "True" ]]; then
                bedtools maskfasta \\
                    -fi <(bowtie2-inspect "{bowtie2_index}" |
                          python "{rename_and_filter_chr}" -f {params.chrom_map} -) \\
                    -bed "{input.mask}" \\
                    -fo >(python "{count_unmasked_bases}" - > "{output}")
            else
                echo {effective_genome_size} > "{output}"
            fi
        }} &> "{log}"
        '''

rule generate_bigwigs:
    input:
        bam = os.path.join(DIR_OUT, "splitbams", "{file}.bam"),
        index = os.path.join(DIR_OUT, "splitbams", "{file}.bam.bai"),
        effective_genome_size = os.path.join(DIR_OUT, "effective_genome_size.txt")
    output:
        os.path.join(DIR_OUT, "bigwigs", "{file}.bw")
    log:
        os.path.join(DIR_LOGS, "bigwig_{file}.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        {{
            if [ ! -s "{input.bam}" ]; then
                # empty BAM file; create empty bigWig file
                echo "BAM file is 0 bytes. Creating empty bigWig file."
                touch "{output}"
            else
                # deepTools bamCoverage currently does not support generating empty bigWig
                # files from BAM files with no aligned reads. See
                # https://github.com/deeptools/deepTools/issues/598
                #
                # This situation can occur when there are clusters with only oligo (BPM) reads
                # and no chromatin (DPM) reads.
                #
                # In such cases, create an empty bigWig file.
                n_reads=$(samtools view -c "{input.bam}")
                if [ $n_reads -eq 0 ]; then
                    echo "- No reads in BAM file for target. Creating empty bigWig file."
                    touch "{output}"
                else

                    if [[ "{bigwig_normalization}" = "RPGC" ]]; then
                        value="$(cat "{input.effective_genome_size}")"
                        effective_genome_size="--effectiveGenomeSize $value"
                    else
                        effective_genome_size=""
                    fi

                    bamCoverage \\
                    --binSize "{binsize}" \\
                    --normalizeUsing "{bigwig_normalization}" \\
                    $effective_genome_size \\
                    -p {threads} \\
                    --bam "{input.bam}" \\
                    --outFileName "{output}"
                fi
            fi
        }} &> "{log}"
        '''
