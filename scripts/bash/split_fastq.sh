#!/bin/bash

####################
# Parse arguments
# $1 number of FASTQ files to split into
# $2 output directory
# $3 output prefix
# $4 number of processes to use for decompression
# $5+ path to input FASTQ files, assumed to be gzipped
####################

if [[ "$#" -lt 5 ]]; then
    echo "Expected 5+ positional arguments."
    echo "Usage: $0 <# chunks> <output directory> <output prefix> <# processes> <input gzipped FASTQ files ...>"
    echo "  Split files are written to <output directory>/<output prefix>##.fastq"
    exit 1
fi

NUM_SPLIT=$1
OUT_DIR=$2
OUT_PREFIX=$3
NUM_PROCESSES=$4
shift 4
INPUT_FASTQS=()
for arg in "$@"; do
    INPUT_FASTQS+=("$arg")
done

####################
# Validate arguments
####################

# check that NUM_SPLIT is a positive integer
if ! [[ "$NUM_SPLIT" =~ ^[1-9][0-9]*$ ]]; then
    echo "The first argument (# of chunks to split into) must be a positive integer."
    exit 1
fi

# check that NUM_PROCESSES is a positive integer
if ! [[ "$NUM_PROCESSES" =~ ^[1-9][0-9]*$ ]]; then
    echo "The fourth argument (# of processes to use) must be a positive integer."
    exit 1
fi

# check that each input FASTQ file exists
for path in "${INPUT_FASTQS[@]}"; do
    if [ ! -f "$path" ]; then
        echo "Input FASTQ file does not exist:" "$path"
        exit 1
    fi
done

####################
# Split FASTQ files
####################

# make output directory if it doesn't exist
mkdir -p "$OUT_DIR"

echo "Splitting the following FASTQ files:" "${INPUT_FASTQS[@]}"

LINES=$(unpigz -c -p "$NUM_PROCESSES" "${INPUT_FASTQS[@]}" | wc -l)
NUM_READS=$((LINES / 4))
echo "Number of reads:" $NUM_READS
READS_PER_CHUNK=$((NUM_READS / NUM_SPLIT + 1)) # +1 because the division floors the value
LINES_PER_CHUNK=$((READS_PER_CHUNK * 4))
echo "Number of reads per chunk:" $READS_PER_CHUNK
unpigz -c -p "$NUM_PROCESSES" "${INPUT_FASTQS[@]}" |
    split -d --additional-suffix='.fastq' -l $LINES_PER_CHUNK - "$OUT_DIR/$OUT_PREFIX"
