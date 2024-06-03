#!/bin/bash
#
# Verify that pipeline output using default configuration parameters on example data
# matches reference output.
#
# Implementation note: due to non-determinism in the pipeline (see README), BAM file
# outputs cannot be compared directly. Consequently, the outputs can only be compared
# at the level of counts of unique read positions on the genome - i.e., upon converting
# BAM files to BED files.

source ~/.bashrc
conda activate chipdip

DIR_OUTPUT="$1"
if [ -z "$DIR_OUTPUT" ]; then
    DIR_OUTPUT="."
fi

# paths to reference BED files
BED_REF_AB1="assets/AB1-A1.bed"
BED_REF_AB2="assets/AB2-A2.bed"

# hashes for merged splitbam output, converted to 3-column BED files, sorted lexicographically
HASH_REF_AB1="ed3ec0eb6c1bdb954921dd5c34efc3a8"
HASH_REF_AB2="d89dbda765c35bdde188dbdde1a1e161"

# check that reference BED files are not corrupted
hash_ab1=$(md5sum "$BED_REF_AB1" | cut -f 1 -d ' ')
[ "$hash_ab1" != "$HASH_REF_AB1" ] && echo "Corrupt reference BED file $BED_REF_AB1" && exit 1
hash_ab2=$(md5sum "$BED_REF_AB2" | cut -f 1 -d ' ')
[ "$hash_ab2" != "$HASH_REF_AB2" ] && echo "Corrupt reference BED file $BED_REF_AB2" && exit 1

# generate BED files from pipeline merged splitbam output
tmpbed1=$(mktemp ./bed1.XXXXX)
tmpbed2=$(mktemp ./bed2.XXXXX)
bedtools bamtobed -i "$DIR_OUTPUT"/workup_test/splitbams/AB1-A1.bam |
    cut -f 1,2,3 |
    sort > "$tmpbed1"
bedtools bamtobed -i "$DIR_OUTPUT"/workup_test/splitbams/AB2-A2.bam |
    cut -f 1,2,3 |
    sort > "$tmpbed2"

diff "$tmpbed1" "$BED_REF_AB1"
# delete BED files if they match reference
[ "$?" = 0 ] && echo "AB1-A1 matches reference" && rm "$tmpbed1"

diff "$tmpbed2" "$BED_REF_AB2"
# delete BED files if they match reference
[ "$?" = 0 ] && echo "AB2-A2 matches reference" && rm "$tmpbed2"
