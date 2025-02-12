#!/bin/bash
set -e

source ~/.bashrc
conda activate chipdip

# assume that this script is run from the working directory set to the tests folder
path_script="../scripts/python/rename_and_filter_chr.py"
path_assets="assets/rename_and_filter_chr"

# chrom maps and inputs tested
#
# chrom_map                                          | input                  | sort option | reference output             
# -------------------------------------------------- | ---------------------- | ----------- | -----------------------------
# re-ordered relative to input: chrom_map_sort.txt   | unsorted: unsorted.sam | auto        | reorder_unsorted_and_sort.sam
# re-ordered relative to input: chrom_map_sort.txt   | unsorted: unsorted.sam | true        | reorder_unsorted_and_sort.sam
# re-ordered relative to input: chrom_map_sort.txt   | unsorted: unsorted.sam | false       | reorder_unsorted.sam         
# same order as input: chrom_map_rename.txt          | unsorted: unsorted.sam | auto        | rename_unsorted.sam          
# same order as input: chrom_map_rename.txt          | unsorted: unsorted.sam | true        | rename_unsorted_and_sort.sam 
# same order as input: chrom_map_rename.txt          | unsorted: unsorted.sam | false       | rename_unsorted.sam          
# re-ordered relative to input: chrom_map_rename.txt | sorted: sorted.sam     | auto        | rename_unsorted_and_sort.sam 
# re-ordered relative to input: chrom_map_rename.txt | sorted: sorted.sam     | true        | rename_unsorted_and_sort.sam 
# re-ordered relative to input: chrom_map_rename.txt | sorted: sorted.sam     | false       | reorder_sorted.sam           
# same order as input: chrom_map_sort.txt            | sorted: sorted2.sam    | auto        | rename_sorted2.sam           
# same order as input: chrom_map_sort.txt            | sorted: sorted2.sam    | true        | reorder_unsorted_and_sort.sam
# same order as input: chrom_map_sort.txt            | sorted: sorted2.sam    | false       | rename_sorted2.sam           
# same order as input: chrom_map_subset.txt          | sorted: sorted.sam     | auto        | rename_subset_sorted.sam     
# same order as input: chrom_map_subset.txt          | sorted: sorted2.sam    | auto        | rename_subset_sorted2.sam    
# same order as input: chrom_map_subset.txt          | sorted: sorted2.sam    | true        | rename_subset_sorted.sam     

# input descriptions
# - sorted.sam: header chromosome names and reads are sorted; header SO tag = "SO:coordinate"
# - sorted2.sam: identical to sorted.sam, except header SO tag = "SO:unknown"
# - unsorted.sam: neither header chromosome names nor reads are sorted

# chrom_map descriptions
# - chrom_map_same.txt: no difference between 2 columns; refX --> refX
# - chrom_map_sort.txt: rename refX to chrX, and chromosome names are sorted
# - chrom_map_rename.txt: rename refX to chrX, but chromosome names are not sorted (1, 3, 2), matching order in unsorted.sam
# - chrom_map_subset.txt: rename refX to chrX, only keep chromosomes chr1 and chr2

# total of 8 unique reference ouputs
# - reorder_unsorted_and_sort.sam, reorder_sorted_and_sort.sam: sorted.sam with chromosomes renamed
# - reorder_unsorted.sam: rename chromosomes and reorder them in the header
# - rename_unsorted_and_sort.sam: renamed chromosomes, and reads sorted to chromosomes in order 1, 3, 2
# - rename_unsorted.sam: renamed chromosomes only
# - reorder_sorted.sam: renamed chromosomes, reordered header, but otherwise identical alignment section to sorted.sam
# - rename_sorted2.sam: identical to reorder_unsorted_and_sort.sam but with header SO tag of "SO:unknown"
# - rename_subset_sorted.sam: rename and filter chromosomes; reads remain sorted
# - rename_subset_sorted2.sam: identical to rename_subset_sorted.sam but with header SO tag of "SO:unknown"

path_unsorted="$path_assets/unsorted.sam"
path_sorted="$path_assets/sorted.sam"
path_sorted2="$path_assets/sorted2.sam"

path_reorder_unsorted_and_sort="$path_assets/reorder_unsorted_and_sort.sam"
path_reorder_unsorted="$path_assets/reorder_unsorted.sam"
path_rename_unsorted_and_sort="$path_assets/rename_unsorted_and_sort.sam"
path_rename_unsorted="$path_assets/rename_unsorted.sam"
path_reorder_sorted="$path_assets/reorder_sorted.sam"
path_rename_sorted2="$path_assets/rename_sorted2.sam"
path_rename_subset_sorted="$path_assets/rename_subset_sorted.sam"
path_rename_subset_sorted2="$path_assets/rename_subset_sorted2.sam"

path_chrom_map_same="$path_assets/chrom_map_same.txt"
path_chrom_map_sort="$path_assets/chrom_map_sort.txt"
path_chrom_map_rename="$path_assets/chrom_map_rename.txt"
path_chrom_map_subset="$path_assets/chrom_map_subset.txt"

path_out="out.bam"

#########
# Tests without using chrom map
#########

path_in="$path_unsorted"

[ -f "$path_out" ] && rm "$path_out"
echo "test 1"
python "$path_script" "$path_in" > "$path_out"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 2"
python "$path_script" -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 3"
python "$path_script" --try-symlink "$path_in" > "$path_out"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 4"
python "$path_script" --try-symlink -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")
[ ! -L "$path_out" ] && echo "no symbolic link at $path_out" && exit 1

#########
# Tests using a chrom map that doesn't perform any renaming or reordering
#########

path_in="$path_sorted"
path_chrom_map="$path_chrom_map_same"

[ -f "$path_out" ] && rm "$path_out"
echo "test 5"
python "$path_script" -c "$path_chrom_map" "$path_in" > "$path_out"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 6"
python "$path_script" -c "$path_chrom_map" -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 7"
python "$path_script" -c "$path_chrom_map" --try-symlink "$path_in" > "$path_out"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")

[ -f "$path_out" ] && rm "$path_out"
echo "test 8"
python "$path_script" -c "$path_chrom_map" --try-symlink -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_in") <(samtools view -h --no-PG "$path_out")
[ ! -L "$path_out" ] && echo "no symbolic link at $path_out" && exit 1

#########
# Tests using subsetting chrom map
#########

path_chrom_map="$path_chrom_map_subset"

path_in="$path_sorted"
path_ref="$path_rename_subset_sorted"

[ -f "$path_out" ] && rm "$path_out"
echo "test 9"
python "$path_script" -c "$path_chrom_map" --try-symlink --sort auto --no-PG -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

path_in="$path_sorted2"
path_ref="$path_rename_subset_sorted2"

[ -f "$path_out" ] && rm "$path_out"
echo "test 10"
python "$path_script" -c "$path_chrom_map" --try-symlink --sort auto --no-PG -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

path_in="$path_sorted2"
path_ref="$path_rename_subset_sorted"

[ -f "$path_out" ] && rm "$path_out"
echo "test 11"
python "$path_script" -c "$path_chrom_map" --try-symlink --sort true --no-PG -o "$path_out" "$path_in"
diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

#########
# Tests using unsorted BAM with re-sorting chrom_map, using default sort option (equivalent to --sort auto)
#########

chrom_maps=(
    "$path_chrom_map_sort"
    "$path_chrom_map_sort"
    "$path_chrom_map_sort"
    "$path_chrom_map_rename"
    "$path_chrom_map_rename"
    "$path_chrom_map_rename"
    "$path_chrom_map_rename"
    "$path_chrom_map_rename"
    "$path_chrom_map_rename"
    "$path_chrom_map_sort"
    "$path_chrom_map_sort"
    "$path_chrom_map_sort"
)

inputs=(
    "$path_unsorted"
    "$path_unsorted"
    "$path_unsorted"
    "$path_unsorted"
    "$path_unsorted"
    "$path_unsorted"
    "$path_sorted"
    "$path_sorted"
    "$path_sorted"
    "$path_sorted2"
    "$path_sorted2"
    "$path_sorted2"
)

references=(
    "$path_reorder_unsorted_and_sort"
    "$path_reorder_unsorted_and_sort"
    "$path_reorder_unsorted"
    "$path_rename_unsorted"
    "$path_rename_unsorted_and_sort"
    "$path_rename_unsorted"
    "$path_rename_unsorted_and_sort"
    "$path_rename_unsorted_and_sort"
    "$path_reorder_sorted"
    "$path_rename_sorted2"
    "$path_reorder_unsorted_and_sort"
    "$path_rename_sorted2"
)

sort_options=( 'auto' 'true' 'false' )

for (( i=0; i<${#chrom_maps[@]}; i++ )); do
    path_chrom_map="${chrom_maps[$i]}"
    path_in="${inputs[i]}"
    path_ref="${references[i]}"

    sort_index=$(( i % 3 ))
    sort_option="${sort_options[$sort_index]}"
    echo "chrom_map: $(basename "$path_chrom_map"), in: $(basename "$path_in"), ref: $(basename "$path_ref"), sort: $sort_option"

    [ -f "$path_out" ] && rm "$path_out"
    echo "- test v1"
    python "$path_script" -c "$path_chrom_map" --sort "$sort_option" -t 4 --no-PG -q "$path_in" > "$path_out" 2> /dev/null
    diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

    [ -f "$path_out" ] && rm "$path_out"
    echo "- test v2"
    python "$path_script" -c "$path_chrom_map" --sort "$sort_option" -t 4 --no-PG -q --try-symlink "$path_in" > "$path_out" 2> /dev/null
    diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

    [ -f "$path_out" ] && rm "$path_out"
    echo "- test v3"
    python "$path_script" -c "$path_chrom_map" --sort "$sort_option" -t 4 --no-PG -q -o "$path_out" "$path_in" 2> /dev/null
    diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")

    [ -f "$path_out" ] && rm "$path_out"
    echo "- test v4"
    python "$path_script" -c "$path_chrom_map" --sort "$sort_option" -t 4 --no-PG -q --try-symlink -o "$path_out" "$path_in" 2> /dev/null
    diff <(samtools view -h --no-PG "$path_ref") <(samtools view -h --no-PG "$path_out")
done
