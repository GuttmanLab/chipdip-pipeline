"""
Aim: A Snakemake workflow to process CHIP-DIP data
"""

import json
import os
import sys
import datetime

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

try:
    email = config["email"]
except:
    email = None
    print("Will not send email on error", file=sys.stderr)

##############################################################################
# Location of scripts
##############################################################################

try:
    DIR_SCRIPTS = config["scripts_dir"]
except:
    print("Scripts directory not specificed in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

split_fastq = os.path.join(DIR_SCRIPTS, "bash/split_fastq.sh")
barcode_id_jar = os.path.join(DIR_SCRIPTS, "java/BarcodeIdentification_v1.2.0.jar")
lig_eff = os.path.join(DIR_SCRIPTS, "python/get_ligation_efficiency.py")
split_bpm_dpm = os.path.join(DIR_SCRIPTS, "python/split_dpm_bpm_fq.py")
validate = os.path.join(DIR_SCRIPTS, "python/validate.py")
rename_and_filter_chr = os.path.join(DIR_SCRIPTS, "python/rename_and_filter_chr.py")
get_clusters = os.path.join(DIR_SCRIPTS, "python/get_clusters.py")
merge_clusters = os.path.join(DIR_SCRIPTS, "python/merge_clusters.py")
fq_to_bam = os.path.join(DIR_SCRIPTS, "python/fastq_to_bam.py")
tag_and_split = os.path.join(DIR_SCRIPTS, "python/threshold_tag_and_split.py")
calculate_effective_genome_size = os.path.join(DIR_SCRIPTS, "python/calculate_effective_genome_size.py")

cluster_statistics = os.path.join(DIR_SCRIPTS, "python/generate_cluster_statistics.py")
cluster_sizes = os.path.join(DIR_SCRIPTS, "python/get_bead_size_distribution.py")
cluster_ecdfs = os.path.join(DIR_SCRIPTS, "python/max_representation_ecdfs_perlib.py")

pipeline_counts = os.path.join(DIR_SCRIPTS, "python/pipeline_counts.py")

##############################################################################
# Load settings
##############################################################################

try:
    bid_config = config["bID"]
    print("Using BarcodeID config: ", bid_config, file=sys.stderr)
except:
    bid_config = "config.txt"
    print('Config "bID" not specified, looking for config at:', bid_config, file=sys.stderr)

formatfile = config.get("format")
if formatfile is not None:
    print("Using split-pool format file:", formatfile, file=sys.stderr)
else:
    print("(WARNING) Format file not specified. The pipeline will NOT ensure barcodes are valid.", file=sys.stderr)

try:
    num_tags = int(config["num_tags"])
    print("Using", num_tags, "tags", file=sys.stderr)
except:
    num_tags = 6
    print('Config "num_tags" not specified, using:', num_tags, file=sys.stderr)

try:
    samples = config["samples"]
    print("Using samples file: ", samples, file=sys.stderr)
except:
    samples = "./samples.json"
    print("Defaulting to working directory for samples JSON file", file=sys.stderr)

try:
    out_dir = config["output_dir"]
    print("All data will be written to: ", out_dir, file=sys.stderr)
except:
    out_dir = os.getcwd()
    print("Defaulting to working directory as output directory", file=sys.stderr)

try:
    temp_dir = config["temp_dir"]
    print("Using temporary directory: ", temp_dir, file=sys.stderr)
except:
    temp_dir = "/central/scratch/"
    print("Defaulting to /central/scratch as temporary directory", file=sys.stderr)

try:
    num_chunks = int(config["num_chunks"])
    print("Splitting FASTQ files into {} chunks for parallel processing".format(num_chunks),
          file=sys.stderr)
except:
    print("Defaulting to 2 chunks for parallel processing", file=sys.stderr)
    num_chunks = 2

try:
    conda_env = config["conda_env"]
except:
    print("No conda environment specified. Defaulting to envs/chipdip.yaml", file=sys.stderr)
    conda_env = "envs/chipdip.yaml"
if conda_env.lower().endswith(".yaml") or conda_env.lower().endswith(".yml"):
    print("Will create new conda environment from", conda_env, file=sys.stderr)
else:
    print("Using existing conda environment:", conda_env, file=sys.stderr)

merge_samples = config.get("merge_samples", False)

##############################################################################
# Load Post Clustering Setting
##############################################################################

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
    print("Will not generate bigWigs as split BAMs are not being generated", file=sys.stderr)
    binsize = False

##############################################################################
# Trimming Sequences
##############################################################################

try:
    adapters = "-g file:" + config["cutadapt_dpm"]
    print("Using cutadapt sequence file", adapters, file=sys.stderr)
except:
    print("DPM adaptor sequences not specificed in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

try:
    oligos = "-g file:" + config["cutadapt_oligos"]
    print("Using bead oligo file", oligos, file=sys.stderr)
except:
    print("Oligo sequences not specified in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

##############################################################################
# DNA Mask
##############################################################################

mask = config.get("mask")
if mask is not None:
    print("Masking reads that align to regions in:", mask, file=sys.stderr)
else:
    mask = ""
    print("(WARNING) Mask path (mask) not specified in config.yaml, no masking will be performed.", file=sys.stderr)

##############################################################################
# Aligner Indexes
##############################################################################

try:
    bowtie2_index = config["bowtie2_index"]
except:
    print("Bowtie2 index not specified in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

##############################################################################
# Chromosomes to keep and/or rename
##############################################################################

path_chrom_map = config.get("path_chrom_map")
if path_chrom_map is None:
    print("Chromosome names not specified, will use all chromosomes in the Bowtie 2 index.",
          file=sys.stderr)

##############################################################################
# Make output directories
##############################################################################

DIR_WORKUP = os.path.join(out_dir, "workup")
DIR_LOGS = os.path.join(DIR_WORKUP, "logs")

DIR_LOGS_CLUSTER = os.path.join(DIR_LOGS, "cluster")
os.makedirs(DIR_LOGS_CLUSTER, exist_ok=True)
out_created = os.path.exists(DIR_LOGS_CLUSTER)
print("Output logs path created:", out_created, file=sys.stderr)

##############################################################################
# Get sample files
##############################################################################

with open(samples) as f:
    FILES = json.load(f)
ALL_SAMPLES = sorted(FILES.keys())

NUM_CHUNKS = [f"{i:03}" for i in range(num_chunks)]

##############################################################################
# Logging
##############################################################################

CONFIG = [os.path.join(DIR_LOGS, "config_" + run_date + ".json")]

LE_LOG_ALL = [os.path.join(DIR_WORKUP, "ligation_efficiency.txt")]

MULTI_QC = [os.path.join(DIR_WORKUP, "qc", "multiqc_report.html")]

LOG_VALIDATE = [os.path.join(DIR_LOGS, "validate.txt")]

##############################################################################
# Trimming
##############################################################################

SPLIT_FQ = expand(
    os.path.join(DIR_WORKUP, "splitfq", "{sample}_{read}.part_{splitid}.fastq.gz"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM = expand(
    [os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz"),
     os.path.join(DIR_WORKUP, "trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

TRIM_LOG = expand(
    os.path.join(DIR_WORKUP, "trimmed/{sample}_{read}.part_{splitid}.fastq.gz_trimming_report.txt"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM_RD = expand(
    [os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz"),
     os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# Barcoding
##############################################################################

BARCODEID = expand(
    os.path.join(DIR_WORKUP, "fastqs/{sample}_{read}.part_{splitid}.barcoded.fastq.gz"),
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

SPLIT_DPM_BPM = expand(
    [os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
     os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz")],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# DNA workup
##############################################################################

Bt2_DNA_ALIGN = expand(
    os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_DNA = expand(
    os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.bam"),
    sample=ALL_SAMPLES)

CHR_DNA = expand(
    os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MASKED = expand(
    os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# Bead workup
##############################################################################

FQ_TO_BAM = expand(
    os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.BPM.bam"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_BEAD = expand(
    os.path.join(DIR_WORKUP, "alignments/{sample}.merged.BPM.bam"),
    sample=ALL_SAMPLES)

##############################################################################
# Clustering
##############################################################################

CLUSTERS = expand(
    os.path.join(DIR_WORKUP, "clusters_parts/{sample}.part_{splitid}.clusters"),
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

CLUSTERS_MERGED = expand(
    os.path.join(DIR_WORKUP, "clusters/{sample}.clusters"),
    sample=ALL_SAMPLES)

##############################################################################
# Post Clustering
##############################################################################

CLUSTER_STATISTICS = [os.path.join(DIR_WORKUP, "clusters/cluster_statistics.txt")]

CLUSTER_SIZES = [os.path.join(DIR_WORKUP, "clusters/DPM_read_distribution.pdf"),
         os.path.join(DIR_WORKUP, "clusters/DPM_cluster_distribution.pdf"),
         os.path.join(DIR_WORKUP, "clusters/BPM_cluster_distribution.pdf"),
         os.path.join(DIR_WORKUP, "clusters/BPM_read_distribution.pdf")]

ECDFS = [os.path.join(DIR_WORKUP, "clusters/Max_representation_ecdf.pdf"),
         os.path.join(DIR_WORKUP, "clusters/Max_representation_counts.pdf")]

SPLITBAMS = expand(
    os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.labeled.bam"),
    sample=ALL_SAMPLES)

SPLITBAMS_STATISTICS = [os.path.join(DIR_WORKUP, "splitbams/splitbam_statistics.txt")]

CLUSTERS_ALL = [os.path.join(DIR_WORKUP, "clusters/all.clusters")]

SPLITBAMS_ALL_LOG = [os.path.join(DIR_LOGS, "splitbams_all.log")]

BIGWIGS_LOG = [os.path.join(DIR_LOGS, "bigwigs.log")]

PIPELINE_COUNTS = [os.path.join(DIR_WORKUP, "pipeline_counts.txt")]

FINAL = \
    MERGE_BEAD + CLUSTER_SIZES + ECDFS + CLUSTER_STATISTICS + MULTI_QC + \
    LE_LOG_ALL + CONFIG + PIPELINE_COUNTS

if binsize:
    FINAL.extend(BIGWIGS_LOG)
else:
    if generate_splitbams and merge_clusters:
        FINAL.extend(SPLITBAMS_ALL_LOG + CLUSTERS_ALL)
    elif generate_splitbams:
        FINAL.extend(SPLITBAMS_STATISTICS)
    elif merge_clusters:
        FINAL.extend(CLUSTERS_ALL + MERGE_DNA)

# ALL_OUTPUTS = \
#     SPLIT_FQ + TRIM + BARCODEID + SPLIT_DPM_BPM + TRIM_RD + \
#     FQ_TO_BAM + MERGE_BEAD + \
#     Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + CLUSTERS + CLUSTERS_MERGED + \
#     CLUSTER_STATISTICS + CLUSTER_SIZES + ECDFS + CLUSTERS_ALL + \
#     SPLITBAMS + SPLITBAMS_STATISTICS + SPLITBAMS_ALL_LOG + BIGWIGS_LOG + \
#     CONFIG + LE_LOG_ALL + MULTI_QC + PIPELINE_COUNTS

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
    shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    sample = "[^\.]+"

# remove all output, leaving just the following in the workup folder:
# - bigwigs/
# - clusters/
# - qc/
# - splitbams/
# - ligation_efficiency.txt
# - pipeline_counts.txt
rule clean:
    shell:
        '''
        for path in {DIR_WORKUP}/*; do
            if [[ "$path" != "{DIR_WORKUP}/bigwigs" ]] &&
               [[ "$path" != "{DIR_WORKUP}/clusters" ]] &&
               [[ "$path" != "{DIR_WORKUP}/qc" ]] &&
               [[ "$path" != "{DIR_WORKUP}/splitbams" ]] &&
               [[ "$path" != "{DIR_WORKUP}/ligation_efficiency.txt" ]] &&
               [[ "$path" != "{DIR_WORKUP}/pipeline_counts.txt" ]]; then
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
            bowtie2-inspect --summary '{bowtie2_index}' > '{log.bt2_sum}'
            python {validate} -c '{input.config}' --bt2_index_summary '{log.bt2_sum}'
        }} &> '{log.log}'
        '''

##############################################################################
# Trimming and barcode identification
##############################################################################

# Split fastq files into chunks to processes in parallel
rule splitfq:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        temp(expand(
            [os.path.join(DIR_WORKUP, "splitfq/{{sample}}_R1.part_{splitid}.fastq"),
             os.path.join(DIR_WORKUP, "splitfq/{{sample}}_R2.part_{splitid}.fastq")],
             splitid=NUM_CHUNKS))
    params:
        dir = os.path.join(DIR_WORKUP, "splitfq"),
        prefix_r1 = "{sample}_R1.part_0",
        prefix_r2 = "{sample}_R2.part_0"
    log:
        os.path.join(DIR_LOGS, "{sample}.splitfq.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        mkdir -p "{params.dir}"
        bash "{split_fastq}" "{input.r1}" {num_chunks} "{params.dir}" "{params.prefix_r1}" {threads} &> {log}
        bash "{split_fastq}" "{input.r2}" {num_chunks} "{params.dir}" "{params.prefix_r2}" {threads} &>> {log}
        '''

# Compress the split fastq files
rule compress_fastq:
    input:
        r1 = os.path.join(DIR_WORKUP, "splitfq/{sample}_R1.part_{splitid}.fastq"),
        r2 = os.path.join(DIR_WORKUP, "splitfq/{sample}_R2.part_{splitid}.fastq")
    output:
        r1 = os.path.join(DIR_WORKUP, "splitfq/{sample}_R1.part_{splitid}.fastq.gz"),
        r2 = os.path.join(DIR_WORKUP, "splitfq/{sample}_R2.part_{splitid}.fastq.gz")
    conda:
        conda_env
    threads:
        8
    shell:
        '''
        pigz -p {threads} "{input.r1}"
        pigz -p {threads} "{input.r2}"
        '''

# Trim adaptors
rule adaptor_trimming_pe:
    input:
        [os.path.join(DIR_WORKUP, "splitfq/{sample}_R1.part_{splitid}.fastq.gz"),
         os.path.join(DIR_WORKUP, "splitfq/{sample}_R2.part_{splitid}.fastq.gz")]
    output:
         os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz"),
         os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.fastq.gz_trimming_report.txt"),
         os.path.join(DIR_WORKUP, "trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"),
         os.path.join(DIR_WORKUP, "trimmed/{sample}_R2.part_{splitid}.fastq.gz_trimming_report.txt")
    params:
        dir = os.path.join(DIR_WORKUP, "trimmed")
    threads:
        10
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.trim_galore.log")
    conda:
        conda_env
    shell:
        '''
        if [[ {threads} -gt 8 ]]; then
            cores=2
        else
            cores=1
        fi

        trim_galore \
        --paired \
        --gzip \
        --cores $cores \
        --quality 20 \
        --fastqc \
        -o "{params.dir}" \
        {input} &> "{log}"
        '''

# Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz"),
        r2 = os.path.join(DIR_WORKUP, "trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz")
    output:
        r1_barcoded = os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz"),
        r2_barcoded = os.path.join(DIR_WORKUP, "fastqs/{sample}_R2.part_{splitid}.barcoded.fastq.gz")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.bID.log")
    shell:
        '''
        java -jar "{barcode_id_jar}" \
        --input1 "{input.r1}" --input2 "{input.r2}" \
        --output1 "{output.r1_barcoded}" --output2 "{output.r2_barcoded}" \
        --config "{bid_config}" &> "{log}"
        '''

# Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        temp(os.path.join(DIR_WORKUP, "{sample}.part_{splitid}.ligation_efficiency.txt"))
    conda:
        conda_env
    shell:
        '''
        python "{lig_eff}" "{input.r1}" "{bid_config}" > "{output}"
        '''

rule cat_ligation_efficiency:
    input:
        expand(
            os.path.join(DIR_WORKUP, "{sample}.part_{splitid}.ligation_efficiency.txt"),
            sample=ALL_SAMPLES,
            splitid=NUM_CHUNKS)
    output:
        LE_LOG_ALL
    shell:
        '''
        tail -n +1 {input} > "{output}"
        '''

# Split barcoded reads into BPM and DPM, remove incomplete barcodes
rule split_bpm_dpm:
    input:
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"),
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_other.fastq.gz"),
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_short.fastq.gz")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.BPM_DPM.log")
    params:
        format = f"--format '{formatfile}'" if formatfile is not None else ""
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
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz")
    output:
        fastq = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz"),
        qc = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.qc.txt")
    params:
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA " + adapters,
        others = "--minimum-length 20"
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.DPM.cutadapt.log")
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         {params.others} \
         -o "{output.fastq}" \
         -j {threads} \
         "{input}" > "{output.qc}") &> "{log}"

        fastqc "{output.fastq}"
        '''

# Trim 9mer oligo sequence from read1 of BPM reads
rule cutadapt_oligo:
    input:
        os.path.join(DIR_WORKUP, "fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz")
    output:
        fastq = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"),
        qc = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt")
    params:
        adapters_r1 = oligos
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.BPM.cutadapt.log")
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         -o "{output.fastq}" \
         -j {threads} \
         "{input}" > "{output.qc}") &> "{log}"
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
        fq = os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz")
    output:
        sorted = os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam"),
        bam = temp(os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.unsorted.bam"))
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.bowtie2.log")
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        (bowtie2 \
         -p 10 \
         -t \
         --phred33 \
         -x "{bowtie2_index}" \
         -U "{input.fq}" | \
         samtools view -bq 20 -F 4 -F 256 - > "{output.bam}") &> "{log}"
        samtools sort -@ {threads} -o "{output.sorted}" "{output.bam}"
        '''

# Rename chromosome names and filter for chromosomes of interest
rule rename_and_filter_chr:
    input:
        os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam"),
    output:
        os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.bam"),
    params:
        chrom_map = f"--chrom_map '{path_chrom_map}'" if path_chrom_map is not None else "",
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.rename_and_filter_chr.log"),
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        python "{rename_and_filter_chr}" {params.chrom_map} -t {threads} "{input}" "{output}" &> "{log}"
        '''

# Merge mask
# - This should increase the speed of the repeat_mask rule compared to a bed file with many overlapping regions.
# - The merged mask is also used in the generate_bigwigs rule, if that is enabled.
rule merge_mask:
    output:
        temp(os.path.join(DIR_WORKUP, "mask_merge.bed"))
    conda:
        conda_env
    shell:
        '''
        if [ -n "{mask}" ]; then
            sort -k1,1 -k2,2n "{mask}" | bedtools merge > "{output}"
        else
            touch "{output}"
        fi
        '''

# Repeat mask aligned DNA reads
rule repeat_mask:
    input:
        bam = os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.bam"),
        mask = os.path.join(DIR_WORKUP, "mask_merge.bed")
    output:
        os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.repeat_mask.log")
    conda:
        conda_env
    shell:
        '''
        {{
            if [ -z "{mask}" ]; then
                echo "No mask file specified, skipping masking."
                ln -s "{input.bam}" "{output}"
            else
                bedtools intersect -v -a "{input.bam}" -b "{input.mask}" > "{output}"
            fi
        }} &> "{log}"
        '''

# Combine all mapped DNA reads into a single bam file per sample
rule merge_dna:
    input:
        expand(
            os.path.join(DIR_WORKUP, "alignments_parts/{{sample}}.part_{splitid}.DNA.chr.masked.bam"),
            splitid=NUM_CHUNKS)
    output:
        os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.bam")
    conda:
        conda_env
    threads:
        10
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_DNA.log")
    shell:
        '''
        samtools merge -@ {threads} "{output}" {input} &> "{log}"
        '''

##############################################################################
# Workup Bead Oligo
##############################################################################

# Convert the BPM FASTQ reads into a BAM file, keeping the UMI
rule fastq_to_bam:
    input:
        os.path.join(DIR_WORKUP, "trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")
    output:
        sorted = os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.BPM.bam"),
        bam = temp(os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.BPM.unsorted.bam"))
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.make_bam.log")
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        python "{fq_to_bam}" --input "{input}" --output "{output.bam}" --config "{bid_config}" &> "{log}"
        samtools sort -@ {threads} -o "{output.sorted}" "{output.bam}"
        '''

# Combine all oligo reads into a single file per sample
rule merge_beads:
    input:
        expand(
            os.path.join(DIR_WORKUP, "alignments_parts/{{sample}}.part_{splitid}.BPM.bam"),
            splitid=NUM_CHUNKS)
    output:
        os.path.join(DIR_WORKUP, "alignments/{sample}.merged.BPM.bam")
    conda:
        conda_env
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_beads.log")
    threads:
        10
    shell:
        '''
        samtools merge -@ {threads} "{output}" {input} &> "{log}"
        '''

##############################################################################
# Make clusters
##############################################################################

# Make clusters from aligned DNA reads and oligo reads
rule make_clusters:
    input:
        dpm = os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam"),
        bpm = os.path.join(DIR_WORKUP, "alignments_parts/{sample}.part_{splitid}.BPM.bam")
    output:
        unsorted = temp(os.path.join(DIR_WORKUP, "clusters_parts/{sample}.part_{splitid}.unsorted.clusters")),
        sorted = os.path.join(DIR_WORKUP, "clusters_parts/{sample}.part_{splitid}.clusters")
    log:
        os.path.join(DIR_LOGS, "{sample}.{splitid}.make_clusters.log")
    threads:
        4
    conda:
        conda_env
    shell:
        '''
        (python "{get_clusters}" \
        -i "{input.bpm}" "{input.dpm}" \
        -o "{output.unsorted}" \
        -n {num_tags}) &> "{log}"

        sort -k 1 -T "{temp_dir}" --parallel={threads} "{output.unsorted}" > "{output.sorted}"
        '''

# Merge clusters from parallel processing into a single cluster file per sample
rule merge_clusters:
    input:
        expand(
            os.path.join(DIR_WORKUP, "clusters_parts/{{sample}}.part_{splitid}.clusters"),
            splitid=NUM_CHUNKS)
    output:
        mega = temp(os.path.join(DIR_WORKUP, "clusters/{sample}.duplicated.clusters")),
        final = os.path.join(DIR_WORKUP, "clusters/{sample}.clusters")
    log:
        os.path.join(DIR_LOGS, "{sample}.merge_clusters.log")
    conda:
       conda_env
    threads:
        10
    shell:
        '''
        sort -k 1 -T "{temp_dir}" --parallel={threads} -m {input} > "{output.mega}"
        python "{merge_clusters}" -i "{output.mega}" -o "{output.final}" &> "{log}"
        '''

##############################################################################
# Profile clusters
##############################################################################

# Generate simple statistics for clusters
rule generate_cluster_statistics:
    input:
        CLUSTERS_ALL + CLUSTERS_MERGED if merge_samples else CLUSTERS_MERGED
    output:
        CLUSTER_STATISTICS
    log:
        os.path.join(DIR_LOGS, "cluster_statistics.log")
    params:
        dir = os.path.join(DIR_WORKUP, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{cluster_statistics}" --directory "{params.dir}" --pattern .clusters \
            > "{output}" 2> "{log}"
        '''

# Generate ecdfs of oligo distribution
rule generate_cluster_ecdfs:
    input:
        CLUSTERS_ALL + CLUSTERS_MERGED if merge_samples else CLUSTERS_MERGED
    output:
        ecdf = os.path.join(DIR_WORKUP, "clusters/Max_representation_ecdf.pdf"),
        counts = os.path.join(DIR_WORKUP, "clusters/Max_representation_counts.pdf")
    log:
        os.path.join(DIR_LOGS, "cluster_ecdfs.log")
    params:
        dir = os.path.join(DIR_WORKUP, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{cluster_ecdfs}" --directory "{params.dir}" --pattern .clusters \
            --xlim 30 &> "{log}"
        '''

# Profile size distribution of clusters
rule get_size_distribution:
    input:
        CLUSTERS_ALL + CLUSTERS_MERGED if merge_samples else CLUSTERS_MERGED
    output:
        dpm = os.path.join(DIR_WORKUP, "clusters/DPM_read_distribution.pdf"),
        dpm2 = os.path.join(DIR_WORKUP, "clusters/DPM_cluster_distribution.pdf"),
        bpm = os.path.join(DIR_WORKUP, "clusters/BPM_read_distribution.pdf"),
        bpm2 = os.path.join(DIR_WORKUP, "clusters/BPM_cluster_distribution.pdf")
    log:
        os.path.join(DIR_LOGS, "size_distribution.log")
    params:
        dir = os.path.join(DIR_WORKUP, "clusters")
    conda:
        conda_env
    shell:
        '''
        python "{cluster_sizes}" --directory "{params.dir}" --pattern .clusters \
            --readtype BPM &> "{log}"
        python "{cluster_sizes}" --directory "{params.dir}" --pattern .clusters \
            --readtype DPM &>> "{log}"
        '''

##############################################################################
# Logging and MultiQC
##############################################################################

# Aggregate metrics using multiqc
rule multiqc:
    input:
        CLUSTERS_MERGED
    output:
        os.path.join(DIR_WORKUP, "qc/multiqc_report.html")
    log:
        os.path.join(DIR_LOGS, "multiqc.log")
    params:
        dir_qc = os.path.join(DIR_WORKUP, "qc")
    conda:
        conda_env
    shell:
        '''
        multiqc -f -o "{params.dir_qc}" "{DIR_WORKUP}" &> "{log}"
        '''

rule pipeline_counts:
    input:
        SPLIT_FQ + TRIM + BARCODEID + SPLIT_DPM_BPM + TRIM_RD + \
        FQ_TO_BAM + MERGE_BEAD + \
        Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + CLUSTERS + \
        CLUSTERS_MERGED + (SPLITBAMS if generate_splitbams else [])
    output:
        csv = os.path.join(DIR_WORKUP, "qc/pipeline_counts.csv"),
        pretty = os.path.join(DIR_WORKUP, "pipeline_counts.txt")
    log:
        os.path.join(DIR_LOGS, "pipeline_counts.log")
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        (python "{pipeline_counts}" \
           --samples "{samples}" \
           -w "{DIR_WORKUP}" \
           -o "{output.csv}" \
           -n {threads} | \
         column -t -s $'\t' > "{output.pretty}") &> "{log}"
        '''

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split:
    input:
        bam = os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.bam"),
        clusters = os.path.join(DIR_WORKUP, "clusters/{sample}.clusters")
    output:
        os.path.join(DIR_WORKUP, "alignments/{sample}.DNA.merged.labeled.bam")
    log:
        os.path.join(DIR_LOGS, "{sample}.splitbams.log")
    params:
        dir_splitbams = os.path.join(DIR_WORKUP, "splitbams")
    conda:
        conda_env
    shell:
        '''
        python "{tag_and_split}" \
         -i "{input.bam}" \
         -c "{input.clusters}" \
         -o "{output}" \
         -d "{params.dir_splitbams}" \
         --min_oligos {min_oligos} \
         --proportion {proportion} \
         --max_size {max_size} \
         --num_tags {num_tags} &> "{log}"
        '''

# Generate summary statistics of individiual bam files
rule generate_splitbam_statistics:
    input:
        SPLITBAMS
    output:
        SPLITBAMS_STATISTICS
    log:
        os.path.join(DIR_LOGS, "splitbam_statistics.log")
    params:
        dir = os.path.join(DIR_WORKUP, "splitbams"),
        samples = [f"'{sample}'" for sample in ALL_SAMPLES]
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            samples=({params.samples})
            for sample in ${{samples[@]}}; do
                for path in "{params.dir}"/"${{sample}}".DNA.merged.labeled*.bam; do
                    count=$(samtools view -@ {threads} -c "$path")
                    echo -e "${{path}}\t${{count}}" >> "{output}"
                done
            done
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
        sort -k 1 -T "{temp_dir}" --parallel={threads} -m {input} > "{output}"
        '''

rule splitbams_all:
    input:
        SPLITBAMS_STATISTICS
    log:
        SPLITBAMS_ALL_LOG
    params:
        dir = os.path.join(DIR_WORKUP, "splitbams"),
        samples = [f"'{sample}'" for sample in ALL_SAMPLES],
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            targets=$(
                ls "{params.dir}"/*.DNA.merged.labeled_*.bam |
                sed -E -e 's/.*\.DNA\.merged\.labeled_(.*)\.bam/\\1/' |
                sort -u
            )
            echo "$targets"
            for target in ${{targets[@]}}; do
                echo "merging BAM files for target $target"
                path_merged_bam="{params.dir}"/"${{target}}.bam"
                samtools merge -f -@ {threads} "$path_merged_bam" "{params.dir}"/*.DNA.merged.labeled_"$target".bam
                samtools index -@ {threads} "$path_merged_bam"

                # Generate summary statistics of merged bam files
                count=$(samtools view -@ {threads} -c "$path_merged_bam")
                echo -e "${{path_merged_bam}}\t${{count}}" >> "{input}"
            done
        }} &> "{log}"
        '''

rule generate_bigwigs:
    input:
        SPLITBAMS_ALL_LOG,
        sorted_merged_mask = os.path.join(DIR_WORKUP, "mask_merge.bed")
    log:
        BIGWIGS_LOG
    params:
        chrom_map = f"--chrom_map '{path_chrom_map}'" if path_chrom_map is not None else "",
        dir_bam = os.path.join(DIR_WORKUP, "splitbams"),
        dir_bigwig = os.path.join(DIR_WORKUP, "bigwigs"),
        binsize = binsize
    conda:
        conda_env
    threads:
        10
    shell:
        '''
        {{
            mkdir -p "{params.dir_bigwig}"

            targets=$(
                ls "{params.dir_bam}"/*.DNA.merged.labeled_*.bam |
                sed -E -e 's/.*\.DNA\.merged\.labeled_(.*)\.bam/\\1/' |
                sort -u
            )
            for target in ${{targets[@]}}; do
                echo "Generating bigWig for target $target"
                path_bam="{params.dir_bam}"/"${{target}}".bam
                path_bigwig="{params.dir_bigwig}"/"${{target}}".bw

                # deepTools bamCoverage currently does not support generating empty bigWig
                # files from BAM files with no aligned reads. See
                # https://github.com/deeptools/deepTools/issues/598
                #
                # This situation can occur when there are clusters with only oligo (BPM) reads
                # and no chromatin (DPM) reads.
                #
                # In such cases, create an empty (i.e., 0 byte) bigWig file.
                n_reads=$(samtools view -c "$path_bam")
                if [ $n_reads -eq "0" ]; then
                    echo "- No reads in BAM file for target. Creating empty bigWig file."
                    touch "$path_bigwig"
                    continue
                fi

                # calculate genome size from header of BAM file
                # subtract regions (from selected chromosomes) in the mask
                effective_genome_size=$(
                    python {calculate_effective_genome_size} "$path_bam" \
                      -m {input.sorted_merged_mask} \
                      {params.chrom_map} \
                      -t {threads}
                )
                echo "- Effective genome size: $effective_genome_size"

                bamCoverage \
                  --binSize {params.binsize} \
                  --normalizeUsing RPGC \
                  --effectiveGenomeSize $effective_genome_size \
                  -p {threads} \
                  --bam "$path_bam" \
                  --outFileName "$path_bigwig"
            done
        }} &> "{log}"
        '''
