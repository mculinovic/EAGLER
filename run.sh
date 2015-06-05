#!/usr/bin/env bash

if [[ $# -ne 4 ]]; then
    echo "usage: $0 <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta> <extensions_output.fasta>"
    exit 1
fi

name="scaffolder"
tmp_dir="./tmp"

if [[ ! -e $tmp_dir ]]; then
    mkdir $tmp_dir
elif [[ ! -d $tmp_dir ]]; then
    echo "$tmp_dir already exists but is not a tmp_directory" 1>&2
fi

time ./debug/$name $1 $2 $3 $4

# rm -r $tmp_dir
