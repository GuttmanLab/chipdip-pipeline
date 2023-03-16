# !/bin/bash

# $1 input fastq
# $2 num split
# $3 output directory
# $4 output prefix

INPUT_FASTQ=$1
NUM_SPLIT=$2
OUT_DIR=$3
OUT_PREFIX=$4

LINES=$( zcat $INPUT_FASTQ | wc -l)
NUM_READS=$(($LINES/4))
echo $NUM_READS
READS_PER_CHUNK=$(($NUM_READS/$NUM_SPLIT))
READS_PER_CHUNK_FUDGE=$(($READS_PER_CHUNK+1))
LINES_PER_CHUNK=$(($READS_PER_CHUNK_FUDGE*4))
echo $READS_PER_CHUNK
zcat $INPUT_FASTQ | split -d --additional-suffix='.fastq' -l $LINES_PER_CHUNK - $OUT_DIR"/"$OUT_PREFIX



