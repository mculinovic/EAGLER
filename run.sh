#!/usr/bin/env bash

if [[ $# -ne 3 ]]; then
    echo "usage: $0 <ont_reads.fasta> <draft_genome.fasta> <output_file.fasta>"
    exit 1
fi

dir="./tmp"

if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

echo "running contig extension..."
time ./debug/main $1 $2 $3

rm -r $dir

echo $line
echo "RESULTS  $3"
