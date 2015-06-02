#!/usr/bin/env bash

# @file run_tests.sh
# @author Marko Culinovic <marko.culinoivc@gmail.com>
# @brief Script used for testing scaffolder methods

# usage
if [[ $# -ne 3 ]]; then
    echo "usage: $0 <ont_reads.fasta> <draft_genome.fasta> <reference_genome.fasta>"
    exit 1
fi

name="scaffolder"
tmp_dir="./tmp"

if [[ ! -e $tmp_dir ]]; then
   mkdir $tmp_dir
elif [[ ! -d $tmp_dir ]]; then
   echo "$tmp_dir already exists but is not a tmp_directory" 1>&2
fi

# Make time only output elapsed time in seconds
TIMEFORMAT=%R


#indexing draft genome
echo "[BWA] indexing $2"
echo "..."
bwa_index="bwa index $2"
# Keep stdout unmolested
exec 3>&1
{ time $bwa_index 1>&3; } 2>&1 | awk '{
    printf "[BWA] indexing finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
exec 3>&-
echo ""

#acquiring number of processor cores for multithreading
num_threads=$(grep -c ^processor /proc/cpuinfo)

#aligning long reads to draft genome
echo "[BWA] aligning reads to draft genome"
echo "..."
bwa_mem="bwa mem -t $num_threads -x pacbio $2 $1"
{ time $bwa_mem > "./tmp/aln.sam"; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


# filenames for scaffolder output
poa_file="./tmp/poa.fasta"
gr_file="./tmp/gr.fasta"
poa_aln_file="./tmp/poa.sam"
gr_aln_file="./tmp/gr.sam"


# running scaffolder with global realign
echo "[$name] extending contigs using global realignment method"
{ time ./release/$name $1 $2 $gr_file; } 2>&1 | awk '{
    printf "[$name] extension finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1


# running scaffolder with poa
echo "[$name] extending contigs using POA consensus method"ž
{ time ./release/$name -p 1 $1 $2 $poa_file; } 2>&1 | awk '{
    printf "[$name] extension finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1

echo "[BWA] indexing $3"
echo "..."
bwa_index="bwa index $3"
# Keep stdout unmolested
exec 3>&1
{ time $bwa_index 1>&3; } 2>&1 | awk '{
    printf "[BWA] indexing finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
exec 3>&-
echo ""


echo "[BWA] aligning global realign extended contigs to reference genome"
echo "..."
bwa_mem="bwa mem -t $num_threads -x pacbio $3 $gr_file"
{ time $bwa_mem > $gr_aln_file; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "[BWA] aligning POA extended contigs to reference genome"
echo "..."
bwa_mem="bwa mem -t $num_threads -x pacbio $3 $poa_file"
{ time $bwa_mem > $poa_aln_file; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "[ANALYSIS] analyzing .sam files"
python3 scripts/extension_analysis.py $3 $poa_aln_file $gr_aln_file
echo "[ANALYSIS] analysis finshed"

rm -rf "./tmp"