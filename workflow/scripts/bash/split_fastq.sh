#!/bin/bash

####################
# Parse arguments
# $1 number of FASTQ files to split into
# $2 output directory
# $3 output prefix
# $4 number of processes to use for decompression
# $5+ path to input FASTQ files, assumed to be gzipped
#
# Dependencies
# - GNU coreutils: csplit, echo, mkdir, split, wc
# - pigz
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

echo "Splitting the following FASTQ files into ${NUM_SPLIT} total chunks:" "${INPUT_FASTQS[@]}"

LINES=$(unpigz -c -p "$NUM_PROCESSES" "${INPUT_FASTQS[@]}" | wc -l)

# check that the total number of lines in the input FASTQ files is a multiple of 4
if ((LINES % 4 > 0)); then
    echo "Error: The total number of lines ($LINES) in the input FASTQ files is not a multiple of 4."
    echo "This may indicate that the input FASTQ files are corrupt."
    exit 1
fi

NUM_READS=$((LINES / 4))
echo "Number of reads:" $NUM_READS

# check that the total number of lines in the input FASTQ files is at least as large as the number of chunks requested
if [ "$NUM_READS" -lt "$NUM_SPLIT" ]; then
    echo "Error: The number of reads is less than the number of chunks requested."
    exit 1
fi

# If the number of reads is not evenly divisible by the number of chunks, increase the
# number of reads per chunk (RPC) by 1. Then check if the number of chunks that would be
# generated using this increased RPC is the same as the number of chunks requested.
READS_PER_CHUNK=$((NUM_READS / NUM_SPLIT))
if ((NUM_READS % NUM_SPLIT > 0)); then
    READS_PER_CHUNK=$((READS_PER_CHUNK + 1))
fi
EFFECTIVE_NUM_SPLIT=$((NUM_READS / READS_PER_CHUNK))
if ((NUM_READS % READS_PER_CHUNK > 0)); then
    EFFECTIVE_NUM_SPLIT=$((EFFECTIVE_NUM_SPLIT + 1))
fi

if [ $EFFECTIVE_NUM_SPLIT -ne $NUM_SPLIT ]; then
    # split --lines cannot create the exact number of chunks requested --> use csplit (which is slower).
    # This can happen when the requesting a large number of chunks for a small number of reads.
    # - Example: 3 chunks requested for 4 reads (16 lines)
    #   - split --lines=4 --> 4 chunks of 1 read
    #   - split --lines=8 --> 2 chunks of 2 reads
    #   - split --lines=12 --> 1 chunk of 3 reads, 1 chunk of 1 read
    echo "Cannot use GNU split to create the exact number of chunks requested. Using slower GNU csplit..."
    READS_PER_CHUNK=$((NUM_READS / NUM_SPLIT))
    echo "Number of reads per chunk:" $READS_PER_CHUNK
    split_points=""
    for ((i = 1; i < NUM_SPLIT; i++)); do
        split_points+=" $((i * READS_PER_CHUNK * 4 + 1))"
    done

    echo "Split points:" $split_points
    unpigz -c -p "$NUM_PROCESSES" "${INPUT_FASTQS[@]}" |
        csplit -f "$OUT_DIR/$OUT_PREFIX" -b "%02d.fastq" - $split_points
else
    echo "Number of reads per chunk:" $READS_PER_CHUNK
    LINES_PER_CHUNK=$((READS_PER_CHUNK * 4))
    unpigz -c -p "$NUM_PROCESSES" "${INPUT_FASTQS[@]}" |
        split --suffix-length=2 -d --additional-suffix='.fastq' -l $LINES_PER_CHUNK - "$OUT_DIR/$OUT_PREFIX"
fi