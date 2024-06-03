#!/bin/bash

# Arguments
# $1 path to input FASTQ, assumed to be gzipped
# $2 number of FASTQ files to split into
# $3 output directory
# $4 output prefix
# $5 number of processes to use for decompression

INPUT_FASTQ=$1
NUM_SPLIT=$2
OUT_DIR=$3
OUT_PREFIX=$4
NUM_PROCESSES=$5

if [[ -z "$NUM_PROCESSES" ]]; then
    NUM_PROCESSES=1
fi

LINES=$(unpigz -c -p "$NUM_PROCESSES" "$INPUT_FASTQ" | wc -l)
NUM_READS=$((LINES / 4))
echo "Number of reads:" $NUM_READS
READS_PER_CHUNK=$((NUM_READS / NUM_SPLIT + 1)) # +1 because the division floors the value
LINES_PER_CHUNK=$((READS_PER_CHUNK * 4))
echo "Number of reads per chunk:" $READS_PER_CHUNK
filter_cmd="pigz -p $NUM_PROCESSES"' > $FILE'
unpigz -c -p "$NUM_PROCESSES" "$INPUT_FASTQ" |
    split -d --additional-suffix='.fastq.gz' -l $LINES_PER_CHUNK --filter="$filter_cmd" - "$OUT_DIR/$OUT_PREFIX"