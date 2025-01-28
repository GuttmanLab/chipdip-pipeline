#!/bin/bash
#
# Verify that pipeline output using default configuration parameters on example data
# matches reference output.
#
# Usage: verify_merged_splitbams_from_example_data.sh [DIR_OUTPUT] [DIR_TEST_ASSETS] [DIR_TEMP]
# - By default, verify_merged_splitbams_from_example_data.sh assumes that it is run from the
#     pipeline directly (i.e., DIR_OUTPUT='.', DIR_TEST_ASSETS='test/assets', DIR_TEMP='.').
# - DIR_OUTPUT: directory containing pipeline output (i.e., the splitbam folder should be
#     located at DIR_OUTPUT/workup/splitbams)
#     (default: .)
# - DIR_TEST_ASSETS: directory containing test assets (e.g., reference BED files)
#     (default: tests/assets)
# - DIR_TEMP: directory for temporary BED files - useful for troubleshooting if the pipeline
#     output does not match the reference output
#     (default: .)
#
# Dependencies: bedtools, tested with version 2.31
#
# Implementation note: due to non-determinism in the pipeline (see README), BAM file
# outputs cannot be compared directly. Consequently, the outputs can only be compared
# at the level of counts of unique read positions on the genome - i.e., upon converting
# BAM files to BED files.

DIR_OUTPUT="$1"
if [ -z "$DIR_OUTPUT" ]; then
    DIR_OUTPUT="."
fi
DIR_TEST_ASSETS="$2"
if [ -z "$DIR_TEST_ASSETS" ]; then
    DIR_TEST_ASSETS="tests/assets"
fi
DIR_TEMP="$3"
if [ -z "$DIR_TEMP" ]; then
    DIR_TEMP="."
fi

# paths to reference BED files
BED_REF_AB1="$DIR_TEST_ASSETS/AB1-A1.bed"
BED_REF_AB2="$DIR_TEST_ASSETS/AB2-A2.bed"

# hashes for merged splitbam output, converted to 3-column BED files, sorted lexicographically
HASH_REF_AB1="ed3ec0eb6c1bdb954921dd5c34efc3a8"
HASH_REF_AB2="d89dbda765c35bdde188dbdde1a1e161"

# hash for expected cluster statistics file
HASH_REF_CLUSTERS="b2efa9824814cca021e287ba36ebb20f"

# check that reference BED files are not corrupted
hash_ab1=$(md5sum "$BED_REF_AB1" | cut -f 1 -d ' ')
[ "$hash_ab1" != "$HASH_REF_AB1" ] && echo "Corrupt reference BED file $BED_REF_AB1" && exit 1
hash_ab2=$(md5sum "$BED_REF_AB2" | cut -f 1 -d ' ')
[ "$hash_ab2" != "$HASH_REF_AB2" ] && echo "Corrupt reference BED file $BED_REF_AB2" && exit 1

# generate BED files from pipeline merged splitbam output
tmpbed1=$(mktemp ./bed1.XXXXX)
tmpbed2=$(mktemp ./bed2.XXXXX)
bedtools bamtobed -i "$DIR_OUTPUT"/workup/splitbams/AB1-A1.bam |
    cut -f 1,2,3 |
    sort > "$tmpbed1"
bedtools bamtobed -i "$DIR_OUTPUT"/workup/splitbams/AB2-A2.bam |
    cut -f 1,2,3 |
    sort > "$tmpbed2"

diff "$tmpbed1" "$BED_REF_AB1"
# delete BED files if they match reference
[ "$?" = 0 ] && echo "AB1-A1 matches reference." && rm "$tmpbed1"

diff "$tmpbed2" "$BED_REF_AB2"
# delete BED files if they match reference
[ "$?" = 0 ] && echo "AB2-A2 matches reference." && rm "$tmpbed2"

# validate cluster file
hash_cluster=$(export LC_ALL=C; sort "$DIR_OUTPUT"/workup/clusters/cluster_statistics.txt | md5sum | cut -f 1 -d ' ')
[ "$hash_cluster" = "$HASH_REF_CLUSTERS" ] && echo "MD5 checksum of cluster_statistics.txt matches reference."
