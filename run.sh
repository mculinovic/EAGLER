#!/usr/bin/env bash

if [[ $# -ne 3 ]]; then
    echo "usage: $0 <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta>"
    exit 1
fi

name="scaffolder"
tmp_dir="./tmp"

if [[ ! -e $tmp_dir ]]; then
    mktmp_dir $tmp_dir
elif [[ ! -d $tmp_dir ]]; then
    echo "$tmp_dir already exists but is not a tmp_directory" 1>&2
fi

echo "running contig extension..."
time ./debug/$name $1 $2 $3

# rm -r $tmp_dir

echo $line
echo "RESULTS  $3"
