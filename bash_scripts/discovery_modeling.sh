#!/usr/bin/env bash


WORKDIR="/home/schwab/repos/kBioReg/SIMULATION"

while read line
do
    mason_genome -l 1000000 -o $WORKDIR/ref.fasta &>/dev/null
    mkdir "$WORKDIR/tmp"
    split_sequence --input $WORKDIR/ref.fasta --parts 1000 --output $WORKDIR/tmp &>/dev/null
    cat $WORKDIR/tmp/*.fa >$WORKDIR/combined.fna
    rm -rf $WORKDIR/tmp
    kbioreg index -k 7 -m na -o tmp_idx -s 644830 $WORKDIR/combined.fna
    kbioreg query -l 1000 tmp_idx.ibf "$line"
done < $1
