# Count all reads from each output file in the pipeline and summarize in a table and tree view.
#
# Implementation overview
# 1. Each output file in the pipeline has a corresponding count file stored in {DIR_OUT}/pipeline_counts/, where the
#    count file:
#    - is a text file containing a single integer value indicating the read count
#    - is named {level}[.{wildcard1}.{wildcard2}...].{ext}.count, where
#      - the wildcards are sorted
#      - {ext} is one of:
#        - fastq-gz: corresponding to a gzip-compressed FASTQ file
#        - bam: corresponding to a BAM file
# 2. The level hierarchy (a directed acyclic graph) is defined in a YAML file
# 3. The Snakefile contains 2 types of rules:
#    i.  rules to count reads from individual files (count_fastq_gz and count_bam)
#    ii. a rule to aggregate all count files into a summary table and a tree view (pipeline_counts)
#          This rule takes as input all expected count files for the entire pipeline.
# 4. The pipeline_counts script parses the pipeline structure, enumerates expected count files, reads the counts,
#    and generates the summary table and tree view.

# This Snakemake module should be included (include:) in the make Snakefile after the following global variables have
# been defined:
# - DIR_OUT: path to output directory
# - pipeline_structure: dict describing the pipeline outputs
# - ALL_SAMPLES: list of sample names
# - TARGETS: list of target names
# - NUM_CHUNKS: list of split IDs
# - FILES: content of samples.json
# - samples: path to samples.json
# - generate_splitbams: boolean indicating whether split BAMs are generated

import itertools
import os
import re
import string

DIR_COUNTS = os.path.join(DIR_OUT, "pipeline_counts")


def validate_pipeline_structure(pipeline: dict):
    """
    Validate that the pipeline structure adheres to the following:
    1. Read output paths end with expected extensions: .fastq.gz, .fq.gz, or .bam
    2. Parent and base levels (if specified) exist in the pipeline structure.

    Args
    - pipeline: mapping from level to a structured information (a dictionary) describing that output. Except for the
        "data" level, the info dict must contain the key "path" mapping to a list of strings that form the path pattern
        for that output.

    Raises
    - ValueError: if any path pattern does not end with a recognized file extension
    """
    for level, info in pipeline.items():
        path = info.get('path')
        if path and not path[-1].endswith(('.fastq.gz', '.fq.gz', '.bam')):
            raise ValueError(f"Invalid file extension in path pattern for level {level}: {path}")
        assert info.get("parent") in (None, *pipeline.keys()), \
            f"Parent {info.get('parent')} of level {level} is not a valid level."
        assert info.get("base") in (None, *pipeline.keys()), \
            f"Parent {info.get('base')} of level {level} is not a valid level."


def path_to_count_ext(path: str) -> str:
    """
    Given a path pattern, return the corresponding count file extension pattern without the leading period.
    """
    if path.endswith(('.fastq.gz', '.fq.gz')):
        return 'fastq-gz.count'
    elif path.endswith('.bam'):
        return 'bam.count'
    else:
        raise ValueError(f"Unknown file extension for path: {path}")


def file_to_count_options(filename: str, pipeline: dict) -> str:
    """
    Given a target count file, return any additional options for counting reads from that file.

    Args
    - filename: count file filename, excluding the .count extension. Must be of the form '{name}[.{wildcard1}.{wildcard2}...]'
    - pipeline: mapping from level to a structured information (a dictionary) describing that output. Except for the
        "data" level, the info dict must contain the key "path" mapping to a list of strings that form the path pattern
        for that output.

    Returns
    - options: additional options for counting reads from the source file
    """
    level = filename.split('.')[0]
    info = pipeline[level]
    return info.get("count_options", "")


def generate_count_files_all(
    pipeline,
    wildcards: dict[str, list[str]],
    DIR_COUNTS: str,
    data_files: dict[str, dict[str, list[str]]] | None = None,
    generate_splitbams: bool = False,
) -> dict[str, list[str]]:
    """
    Generate all expected count file paths for the entire pipeline.

    Args
    - pipeline: mapping from level to a structured information (a dictionary) describing that output. Except for the
        "data" level, the info dict must contain the key "path" mapping to a list of strings that form the path pattern
        for that output.
    - wildcards: mapping of placeholder names to lists of values
    - data_files: mapping from sample names to a list of their source data files. If None, defaults to the global FILES.
    - DIR_COUNTS: directory to store count files
    - generate_splitbams: whether split BAM files are generated in the pipeline

    Returns
    - path_counts_all: mapping from level to a list of expected count file paths
    """
    if data_files is None:
        data_files = FILES
    formatter = string.Formatter()
    assert sorted(wildcards['sample']) == sorted(data_files.keys())
    path_counts_all = dict()
    for level, info in pipeline.items():
        path_counts = []
        if level == "data":
            for sample, d in data_files.items():
                files_R1 = d["R1"]
                for file_number in range(len(files_R1)):
                    path_counts.append(os.path.join(DIR_COUNTS, f'data.{sample}.{file_number}.fastq-gz.count'))
            path_counts_all[level] = path_counts
            continue
        if level.startswith("splitbams") and not generate_splitbams:
            continue
        wildcards_filtered = dict()
        for field, values in wildcards.items():
            exclude = set(info.get("exclude", dict()).get(field, []))
            wildcards_filtered[field] = [v for v in values if v not in exclude]
        fields = sorted(field for _, field, _, _ in formatter.parse(os.path.join(*info['path'])) if field is not None)
        assert all(field in wildcards_filtered.keys() for field in fields), \
            f"Not all fields {fields} for output level {level} are in wildcards (keys: {wildcards_filtered.keys()})."
        path_count_template = os.path.join(
            DIR_COUNTS,
            '.'.join((
                level,
                '.'.join([f"{{{field}}}" for field in fields]),
                path_to_count_ext(info["path"][-1])
            ))
        )
        for field_combination in itertools.product(*[wildcards_filtered[field] for field in fields]):
            path_counts.append(path_count_template.format(**dict(zip(fields, field_combination))))
        path_counts_all[level] = path_counts
    return path_counts_all

def count_filename_to_source(
    filename: str,
    pipeline: dict,
    source_prefix: str = "",
    data_files: dict[str, dict[str, list[str]]] | None = None,
) -> str:
    """
    Given a target count file, determine the reads files from which to generate the count file.

    Args
    - filename: count file filename, excluding the .count extension. Must be of the form '{name}[.{wildcard1}.{wildcard2}...]'
    - pipeline: mapping from level to a structured information (a dictionary) describing that output. Except for the
        "data" level, the info dict must contain the key "path" mapping to a list of strings that form the path pattern
        for that output.
    - data_files: mapping from sample names to a list of their source data files. If None, defaults to the global FILES.
    - source_prefix: prefix to add to the source path (e.g., output directory)

    Returns
    - path_source: paths to the source reads files
    """
    if data_files is None:
        data_files = FILES

    parts = filename.split('.')
    level = parts[0]
    field_values = parts[1:]
    info = pipeline[level]

    # special case
    if level == "data":
        sample, file_number = field_values
        return data_files[sample]["R1"][int(file_number)]
            
    formatter = string.Formatter()
    fields = sorted(field for _, field, _, _ in formatter.parse(os.path.join(*info['path'])) if field is not None)
    assert len(fields) == len(field_values), f"Filename {filename} does not match expected fields for output level {level}."
    path_source = os.path.join(source_prefix, os.path.join(*info['path']).format(**dict(zip(fields, field_values))))
    return path_source


validate_pipeline_structure(pipeline_structure)
COUNT_FILES_ALL = sorted(set(itertools.chain.from_iterable(generate_count_files_all(
    pipeline_structure,
    wildcards=dict(sample=ALL_SAMPLES, target=TARGETS, splitid=NUM_CHUNKS),
    data_files=FILES,
    DIR_COUNTS=os.path.join(DIR_OUT, "pipeline_counts"),
    generate_splitbams=generate_splitbams
).values())))

# Count reads from a list of gzip-compressed FASTQ files.
rule count_fastq_gz:
    input:
        lambda w: count_filename_to_source(w.file, pipeline_structure, source_prefix=DIR_OUT)
    output:
        temp(os.path.join(DIR_OUT, "pipeline_counts", "{file}.fastq-gz.count"))
    conda:
        conda_env
    threads:
        2
    shell:
        """
        total_lines=$(unpigz -c "{input}" | wc -l)
        remainder=$(( total_lines % 4 ))
        if [ "$remainder" -ne 0 ]; then
            echo "ERROR: FASTQ file "{input}" has $total_lines lines, which is not divisible by 4."
            echo "This suggests the file may be truncated or corrupted."
            exit 1
        fi
        read_count=$(( total_lines / 4 ))
        echo "$read_count" > "{output}"
        """

rule count_bam:
    input:
        lambda w: count_filename_to_source(w.file, pipeline_structure, source_prefix=DIR_OUT)
    output:
        temp(os.path.join(DIR_OUT, "pipeline_counts", "{file}.bam.count"))
    conda:
        conda_env
    params:
        options = lambda w: file_to_count_options(w.file, pipeline_structure)
    shell:
        """
        samtools view -c {params.options} "{input}" > "{output}"
        """

rule pipeline_counts:
    input:
        COUNT_FILES_ALL
    output:
        csv = os.path.join(DIR_OUT, "qc", "pipeline_counts.csv"),
        pretty = os.path.join(DIR_OUT, "pipeline_counts.txt")
    log:
        os.path.join(DIR_LOGS, "pipeline_counts.log")
    conda:
        conda_env
    params:
        samples = samples,
        splitids = NUM_CHUNKS,
        targets = TARGETS,
        dir_counts = DIR_COUNTS,
        path_pipeline = path_pipeline_structure,
    shell:
        '''
        {{
            python {pipeline_counts} \\
                --path_samples {params.samples:q} \\
                --splitids {params.splitids:q} \\
                --targets {params.targets:q} \\
                --dir_counts "{params.dir_counts}" \\
                --out "{output.csv}" \\
                "{params.path_pipeline}" |
            column -t -s $'\\t' > "{output.pretty}"
        }} &> "{log}"
        '''
