"""
Aim: A Snakemake workflow to process CHIP-DIP data
"""

import os
import sys
import numpy as np
import datetime

##############################################################################
# Initialize settings
##############################################################################

# Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime("%Y.%m.%d.")

try:
    config_path = config["config_path"]
except:
    config_path = "config.yaml"

configfile: config_path

try:
    email = config["email"]
except:
    email = None
    print("Will not send email on error", file=sys.stderr)

##############################################################################
# Location of scripts
##############################################################################

barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_bpm_dpm = "scripts/python/split_dpm_bpm_fq.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
merge_clusters = "scripts/python/merge_clusters.py"
fq_to_bam = "scripts/python/fastq_to_bam.py"
split_fastq = "scripts/bash/split_fastq.sh"

cluster_counts = "scripts/python/generate_cluster_statistics.py"
cluster_sizes = "scripts/python/get_bead_size_distribution.py"
cluster_ecdfs = "scripts/python/max_representation_ecdfs_perlib.py"

tag_and_split = "scripts/python/threshold_tag_and_split.py"

##############################################################################
# Load settings
##############################################################################

try:
    bid_config = config["bID"]
    print("Using BarcodeID config: ", bid_config, file=sys.stderr)
except:
    bid_config = "config.txt"
    print('Config "bID" not specified, looking for config at:', bid_config, file=sys.stderr)

try:
    formatfile = config["format"]
    print("Using split-pool format file: ", formatfile, file=sys.stderr)
except:
    formatfile = "format.txt"
    print("Format file not specified, looking for file at:", formatfile, file=sys.stderr)

try:
    num_tags = int(config["num_tags"])
    print("Using", num_tags, "tags", file=sys.stderr)
except:
    num_tags = 6
    print('Config "num_tags" not specified, using:', num_tags, file=sys.stderr)

try:
    assembly = config["assembly"]
    assert assembly in ["mm10", "hg38"], 'Only "mm10" or "hg38" currently supported'
    print("Using", assembly, file=sys.stderr)
except:
    print('Config "assembly" not specified, defaulting to "mm10"', file=sys.stderr)
    assembly = "mm10"

try:
    samples = config["samples"]
    print("Using samples file: ", samples, file=sys.stderr)
except:
    samples = "./samples.json"
    print("Defaulting to working directory for samples json file", file=sys.stderr)

try:
    out_dir = config["output_dir"]
    print("All data will be written to: ", out_dir, file=sys.stderr)
except:
    out_dir = ""
    print("Defaulting to working directory as output directory", file=sys.stderr)

try:
    temp_dir = config["temp_dir"]
    print("Using temporary directory: ", temp_dir, file=sys.stderr)
except:
    temp_dir = "/central/scratch/"
    print("Defaulting to central scratch as temporary directory", file=sys.stderr)

try:
    num_chunks = int(config["num_chunks"])
except:
    num_chunks = 2

##############################################################################
# Load Post Clustering Setting
##############################################################################

try:
    generate_splitbams = config["generate_splitbams"]
except:
    generate_splitbams = False

try:
    min_oligos = config["min_oligos"]
except:
    min_oligos = 2

try:
    proportion = config["proportion"]
except:
    proportion = 0.8

try:
    max_size = config["max_size"]
except:
    max_size = 10000

if generate_splitbams:
    print("Will generate bam files for individual targets using:", file=sys.stderr)
    print("\t min_oligos: ", min_oligos, file=sys.stderr)
    print("\t proportion: ", proportion, file=sys.stderr)
    print("\t max_size: ", max_size, file=sys.stderr)
else:
    print("Will not generate bam files for individual targets.", file=sys.stderr)

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

try:
    mask = config["mask"][config["assembly"]]
except:
    print("Mask path not specified in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

##############################################################################
# Aligner Indexes
##############################################################################

try:
    bowtie2_index = config["bowtie2_index"][config["assembly"]]
except:
    print("Bowtie2 index not specified in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

##############################################################################
# Make output directories
##############################################################################

os.makedirs(out_dir + "workup/logs/cluster", exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
print("Output logs path created:", out_created, file=sys.stderr)

##############################################################################
# Get sample files
##############################################################################

FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get("R1")])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get("R2")])

NUM_CHUNKS = [f"{i:03}" for i in np.arange(0, num_chunks)]

##############################################################################
# Logging
##############################################################################

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]

MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]

##############################################################################
# Trimming
##############################################################################

SPLIT_FQ = expand(
    out_dir + "workup/splitfq/{sample}_{read}.part_{splitid}.fastq.gz",
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM = expand(
    [out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
     out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

TRIM_LOG = expand(
    out_dir + "workup/trimmed/{sample}_{read}.part_{splitid}.fastq.gz_trimming_report.txt",
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

TRIM_RD = expand(
    [out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz",
     out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# Barcoding
##############################################################################

BARCODEID = expand(
    out_dir + "workup/fastqs/{sample}_{read}.part_{splitid}.barcoded.fastq.gz",
    sample=ALL_SAMPLES,
    read=["R1", "R2"],
    splitid=NUM_CHUNKS)

SPLIT_DPM_BPM = expand(
    [out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
     out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"],
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# DNA workup
##############################################################################

Bt2_DNA_ALIGN = expand(
    out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam",
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_DNA = expand(
    out_dir + "workup/alignments/{sample}.DNA.merged.bam",
    sample=ALL_SAMPLES)

CHR_DNA = expand(
    out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.bam",
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MASKED = expand(
    out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam",
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

##############################################################################
# Bead workup
##############################################################################

FQ_TO_BAM = expand(
    out_dir + "workup/alignments_parts/{sample}.part_{splitid}.BPM.bam",
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

MERGE_BEAD = expand(
    out_dir + "workup/alignments/{sample}.merged.BPM.bam",
    sample=ALL_SAMPLES)

##############################################################################
# Clustering
##############################################################################

CLUSTERS = expand(
    out_dir + "workup/clusters_parts/{sample}.part_{splitid}.clusters",
    sample=ALL_SAMPLES,
    splitid=NUM_CHUNKS)

CLUSTERS_MERGED = expand(
    out_dir + "workup/clusters/{sample}.clusters",
    sample=ALL_SAMPLES)

##############################################################################
# Post Clustering
##############################################################################

COUNTS = [out_dir + "workup/clusters/cluster_statistics.txt"]

SIZES = [out_dir + "workup/clusters/DPM_read_distribution.pdf",
         out_dir + "workup/clusters/DPM_cluster_distribution.pdf",
         out_dir + "workup/clusters/BPM_cluster_distribution.pdf",
         out_dir + "workup/clusters/BPM_read_distribution.pdf"]

ECDFS = [out_dir + "workup/clusters/Max_representation_ecdf.pdf",
         out_dir + "workup/clusters/Max_representation_counts.pdf"]

SPLITBAMS = expand(
    out_dir + "workup/alignments/{sample}.DNA.merged.labeled.bam",
    sample=ALL_SAMPLES)

SPLITBAMS_COUNTS = [out_dir + "workup/splitbams/splitbam_statistics.txt"]

if generate_splitbams == False:
    SPLITBAMS = []
    SPLITBAMS_COUNTS = []

##############################################################################
##############################################################################
# RULE ALL
##############################################################################
##############################################################################

rule all:
    input: CONFIG + SPLIT_FQ + ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + SPLIT_DPM_BPM +  MERGE_BEAD + FQ_TO_BAM + Bt2_DNA_ALIGN + CHR_DNA + MASKED + MERGE_DNA + CLUSTERS + CLUSTERS_MERGED + MULTI_QC + COUNTS + SIZES + ECDFS + SPLITBAMS + SPLITBAMS_COUNTS

# Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    sample = "[^\.]+"

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
            [out_dir + "workup/splitfq/{{sample}}_R1.part_{splitid}.fastq",
             out_dir + "workup/splitfq/{{sample}}_R2.part_{splitid}.fastq"],
             splitid=NUM_CHUNKS))
    params:
        dir = out_dir + "workup/splitfq",
        prefix_r1 = "{sample}_R1.part_0",
        prefix_r2 = "{sample}_R2.part_0"
    log:
        out_dir + "workup/logs/{sample}.splitfq.log"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    shell:
        '''
        mkdir -p {params.dir}
        bash {split_fastq} {input.r1} {num_chunks} {params.dir} {params.prefix_r1}
        bash {split_fastq} {input.r2} {num_chunks} {params.dir} {params.prefix_r2}
        '''

# Compress the split fastq files
rule compress_fastq:
    input:
        r1 = out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq",
        r2 = out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq"
    output:
        r1 = out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq.gz",
        r2 = out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq.gz"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    shell:
        '''
        pigz -p {threads} {input.r1}
        pigz -p {threads} {input.r2}
        '''

# Trim adaptors
rule adaptor_trimming_pe:
    input:
        [out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq.gz",
         out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq.gz"]
    output:
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.fastq.gz_trimming_report.txt"
    threads:
        10
    log:
        out_dir + "workup/logs/{sample}.{splitid}.trim_galore.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        if [[ {threads} -gt 8 ]]
        then
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
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''

# Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"
    output:
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"

# Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz"
    output:
        temp(out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt")
    conda:
        "envs/sprite.yaml"
    shell:
        "python {lig_eff} {input.r1} > {output}"

rule cat_ligation_efficiency:
    input:
        expand(
            out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt",
            sample=ALL_SAMPLES,
            splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"

# Split barcoded reads into BPM and DPM, remove incomplete barcodes
rule split_bpm_dpm:
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_other.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.BPM_DPM.log"
    conda:
       "envs/sprite.yaml"
    shell:
        "python {split_bpm_dpm} --r1 {input} --format {formatfile} &> {log}"

##############################################################################
# Cutadapt
##############################################################################

# Trim DPM from read1 of DPM reads, remove DPM dimer reads
rule cutadapt_dpm:
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"
    output:
        fastq = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz",
        qc = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.qc.txt"
    params:
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA " + adapters,
        others = "--minimum-length 20"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.DPM.cutadapt.log"
    threads: 10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         {params.others} \
         -o {output.fastq} \
         -j {threads} \
         {input} > {output.qc}) &> {log}

        fastqc {output.fastq}
        '''

# Trim 9mer oligo sequence from read1 of BPM reads
rule cutadapt_oligo:
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_bpm.fastq.gz"
    output:
        fastq = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz",
        qc = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt"
    params:
        adapters_r1 = oligos
    log:
        out_dir + "workup/logs/{sample}.{splitid}.BPM.cutadapt.log"
    threads: 10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         -o {output.fastq} \
         -j {threads} \
         {input} > {output.qc}) &> {log}
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
        fq = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_dpm.RDtrim.fastq.gz"
    output:
        sorted = out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam",
        bam = temp(out_dir + "workup/alignments_parts/{sample}.part_{splitid}.unsorted.bam")
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bowtie2.log"
    threads:
        10
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (bowtie2 \
         -p 10 \
         -t \
         --phred33 \
         -x {bowtie2_index} \
         -U {input.fq} | \
         samtools view -bq 20 -F 4 -F 256 - > {output.bam}) &> {log}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        '''

# Add 'chr' to chromosome names
rule add_chr:
    input:
        out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam",
    output:
        out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.bam",
    log:
        out_dir + "workup/logs/{sample}.{splitid}.add_chr.log",
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {assembly} &> {log}
        '''

# Repeat mask aligned DNA reads
rule repeat_mask:
    input:
        out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.bam"
    output:
        out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''

# Combine all mapped DNA reads into a single bam file per sample
rule merge_dna:
    input:
        expand(
            out_dir + "workup/alignments_parts/{{sample}}.part_{splitid}.DNA.chr.masked.bam",
            splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{sample}.DNA.merged.bam"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    log:
        out_dir + "workup/logs/{sample}.merge_DNA.log"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) &> {log}
        '''

##############################################################################
# Workup Bead Oligo
##############################################################################

# Convert the BPM FASTQ reads into a BAM file, keeping the UMI
rule fastq_to_bam:
    input:
        out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"
    output:
        sorted = out_dir + "workup/alignments_parts/{sample}.part_{splitid}.BPM.bam",
        bam = temp(out_dir + "workup/alignments_parts/{sample}.part_{splitid}.BPM.unsorted.bam")
    log:
        out_dir + "workup/logs/{sample}.{splitid}.make_bam.log"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    shell:
        '''
        python {fq_to_bam} --input {input} --output {output.bam} --config {bid_config} &> {log}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        '''

# Combine all oligo reads into a single file per sample
rule merge_beads:
    input:
        expand(
            out_dir + "workup/alignments_parts/{{sample}}.part_{splitid}.BPM.bam",
            splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{sample}.merged.BPM.bam"
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{sample}.merge_beads.log"
    threads:
        8
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) >& {log}
        '''

##############################################################################
# Make clusters
##############################################################################

# Make clusters from aligned DNA reads and oligo reads
rule make_clusters:
    input:
        dpm = out_dir + "workup/alignments_parts/{sample}.part_{splitid}.DNA.chr.masked.bam",
        bpm = out_dir + "workup/alignments_parts/{sample}.part_{splitid}.BPM.bam"
    output:
        unsorted = temp(out_dir + "workup/clusters_parts/{sample}.part_{splitid}.unsorted.clusters"),
        sorted = out_dir + "workup/clusters_parts/{sample}.part_{splitid}.clusters"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.make_clusters.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (python {get_clusters} \
        -i {input.bpm} {input.dpm}\
        -o {output.unsorted} \
        -n {num_tags})  &> {log}

        sort -k 1 -T {temp_dir} {output.unsorted} > {output.sorted}
        '''

# Merge clusters from parallel processing into a single cluster file per sample
rule merge_clusters:
    input:
        expand(
            out_dir + "workup/clusters_parts/{{sample}}.part_{splitid}.clusters",
            splitid=NUM_CHUNKS)
    output:
        mega = temp(out_dir + "workup/clusters/{sample}.duplicated.clusters"),
        final = out_dir + "workup/clusters/{sample}.clusters"
    conda:
       "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{sample}.merge_clusters.log"
    shell:
        '''
        sort -k 1 -T {temp_dir} -m {input} > {output.mega}
        (python {merge_clusters} -i {output.mega} -o {output.final}) &> {log}
        '''

##############################################################################
# Profile clusters
##############################################################################

# Generate simple statistics for clusters
rule generate_cluster_statistics:
    input:
        expand([out_dir + "workup/clusters/{sample}.clusters"], sample=ALL_SAMPLES)
    output:
        out_dir + "workup/clusters/cluster_statistics.txt"
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {cluster_counts} --directory {params.dir} --pattern .clusters > {output}
        '''

# Generate ecdfs of oligo distribution
rule generate_cluster_ecdfs:
    input:
        expand([out_dir + "workup/clusters/{sample}.clusters"], sample=ALL_SAMPLES)
    output:
        ecdf = out_dir + "workup/clusters/Max_representation_ecdf.pdf",
        counts = out_dir + "workup/clusters/Max_representation_counts.pdf"
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/plotting.yaml"
    shell:
        '''
        python {cluster_ecdfs} --directory {params.dir} --pattern .clusters --xlim 30
        '''

# Profile size distribution of clusters
rule get_size_distribution:
    input:
        expand([out_dir + "workup/clusters/{sample}.clusters"], sample=ALL_SAMPLES)
    output:
        dpm = out_dir + "workup/clusters/DPM_read_distribution.pdf",
        dpm2 = out_dir + "workup/clusters/DPM_cluster_distribution.pdf",
        bpm = out_dir + "workup/clusters/BPM_read_distribution.pdf",
        bpm2 = out_dir + "workup/clusters/BPM_cluster_distribution.pdf"
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {cluster_sizes} --directory {params.dir} --pattern .clusters --readtype BPM
        python {cluster_sizes} --directory {params.dir} --pattern .clusters --readtype DPM
        '''

##############################################################################
# Logging and MultiQC
##############################################################################

# Copy config.yaml into logs folder with run date
rule log_config:
    input:
        config_path
    output:
        out_dir + "workup/logs/config_" + run_date + "yaml"
    shell:
        "cp {input} {output}"

# Aggregate metrics using multiqc
rule multiqc:
    input:
        expand([out_dir + "workup/clusters/{sample}.clusters"], sample=ALL_SAMPLES)
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda:
        "envs/sprite.yaml"
    shell:
        "multiqc {out_dir}workup -o {out_dir}workup/qc"

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split:
    input:
        bam = out_dir + "workup/alignments/{sample}.DNA.merged.bam",
        clusters = out_dir + "workup/clusters/{sample}.clusters"
    output:
        bam = out_dir + "workup/alignments/{sample}.DNA.merged.labeled.bam",
        touch = temp(touch(out_dir + "workup/splitbams/{sample}.done"))
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{sample}.splitbams.log"
    shell:
        '''
        (python {tag_and_split} \
         -i {input.bam} \
         -c {input.clusters} \
         -o {output.bam} \
         -d workup/splitbams \
         --min_oligos {min_oligos} \
         --proportion {proportion} \
         --max_size {max_size} \
         --num_tags {num_tags}) &> {log}
        '''

# Generate summary statistics of individiual bam files
rule generate_splitbam_statistics:
    input:
        expand([out_dir + "workup/splitbams/{sample}.done"], sample=ALL_SAMPLES)
    output:
        out_dir + "workup/splitbams/splitbam_statistics.txt"
    params:
        dir = out_dir + "workup/splitbams"
    conda:
        "envs/sprite.yaml"
    shell:
        "for f in {params.dir}/*bam; do echo $f; samtools view -c $f; done > {output}"
