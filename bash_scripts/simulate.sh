#!/usr/bin/env bash

# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

set -Eeuxo pipefail

BINARY_DIR="/home/vincent/Desktop//kBioReg/data/bin" # Dir containing "mason_genome", "split_sequence", etc.
OUT_DIR="/home/vincent/Desktop//kBioReg/build" # Where simulated data should be stored
LENGTH=$1 #419430400 # 4*2^30 =  4GiB
SEED=42 # was 20181406 before, but was hardcoded to 42 in seqan
BIN_NUMBER=$2 #64
ERRORS=2
READ_LENGTHS="100 150 250"
READ_COUNT=1048576
HAPLOTYPE_COUNT=1

output_dir=$OUT_DIR/$BIN_NUMBER
bin_dir=$output_dir/bins
info_dir=$output_dir/info

mkdir -p $output_dir
mkdir -p $bin_dir
mkdir -p $info_dir

bin_length=$((LENGTH / BIN_NUMBER))
echo "Simulating $BIN_NUMBER bins with reference length of $LENGTH and bin_length of $bin_length">&2

# Simulate reference
echo "Simulating genome">&2
$BINARY_DIR/mason_genome -l $LENGTH -o $bin_dir/ref.fasta -s $SEED &>/dev/null
# Evenly distribute it over bins
echo "Splitting genome into bins">&2
$BINARY_DIR/split_sequence --input $bin_dir/ref.fasta --length $bin_length --parts $BIN_NUMBER
# We do not need the reference anymore
rm $bin_dir/ref.fasta
#add words to fasta file
MAXCOUNT=99
count=0
num=$(($2-1))
declare -a number=()
while [ "$count" -le $MAXCOUNT ]; do
    number[count]=$(($RANDOM % num ))
    let "count += 1"
done

echo "${number[*]}"
sorted=($(printf '%s\n' "${number[@]}"|sort -n))

sorted[100]=$2

echo "${sorted[*]}"
echo "add words">&2
> bins.txt
occ=0
count=0
index=0
for i in $bin_dir/*.fa
do
    #if [ $index -le 99 ]
     #   then
        while [ $count == ${sorted[$index]} ]; do
            occ=$((occ+1))
            HITWORD=$(shuf -n 1 $OUT_DIR/words.txt)
            echo $HITWORD >> $i
            echo "$i,$HITWORD" >> bins.txt
            let "index += 1"
        done
    #fi
    let "count += 1"
done

# echo "add words">&2
# > bins.txt
# occ=0
# all=0
# for i in $bin_dir/*.fa
# do
#     rand=$((1 + $RANDOM % 100))
#     if [ $rand -le $3 ]
#         then
#             occ=$((occ+1))
#             HITWORD=$(shuf -n 1 $OUT_DIR/words.txt)
#             echo $HITWORD >> $i
#             echo "$i,$HITWORD" >> bins.txt
#     fi
#     all=$((all+1))
# done
echo "Occurence: $occ">&2
echo "all: $count">&2
echo "rate: $(($occ*100 / $count))">&2
cat $bin_dir/*.fa>$bin_dir/all_bins.fa
# rm -r 64/ words.txt bins.txt hits.txt kbioreg_output.log benchmark_idx.ibf benchmark.txt