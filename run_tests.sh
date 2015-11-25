#!/usr/bin/env bash

# @file run_tests.sh
# @author Marko Culinovic <marko.culinoivc@gmail.com>
# @brief Script used for testing scaffolder methods

# usage
if [[ $# -ne 4 ]]; then
    echo "usage: $0 <long_reads.fasta> <draft_genome.fasta> <reference_genome.fasta> <output_dir>"
    exit 1
fi

name="eagler"
tmp_dir="./tmp"

if [[ ! -e $tmp_dir ]]; then
   mkdir $tmp_dir
elif [[ ! -d $tmp_dir ]]; then
   echo "$tmp_dir already exists but is not a tmp_directory" 1>&2
fi

output_dir=$4
if [[ ! -e $output_dir ]]; then
    mkdir $output_dir
fi

bwa_dir="$output_dir/bwa/"
if [[ ! -e $bwa_dir ]]; then
    mkdir $bwa_dir
fi

graphmap_dir="$output_dir/graphmap/"
if [[ ! -e $graphmap_dir ]]; then
    mkdir $graphmap_dir
fi

poa_dir="$output_dir/poa/"
if [[ ! -e $poa_dir ]]; then
    mkdir $poa_dir
fi


# Make time only output elapsed time in seconds
TIMEFORMAT=%R


#indexing draft genome
#echo "[BWA] indexing $2"
#echo "..."
#bwa_index="bwa index $2"
# Keep stdout unmolested
#exec 3>&1
# time $bwa_index 1>&3; } 2>&1 | awk '{
#    printf "[BWA] indexing finished in %d hours, %d minutes and %.3f seconds\n",
#           $1/3600, $1%3600/60, $1%60
#}' | tail -1
#exec 3>&-
#echo ""

#acquiring number of processor cores for multithreading
num_threads=$(grep -c ^processor /proc/cpuinfo)

#aligning long reads to draft genome
#echo "[BWA] aligning reads to draft genome"
#echo "..."
#bwa_mem="bwa mem -t $num_threads -x pacbio -Y $2 $1"
# { time $bwa_mem > "./tmp/aln.sam"; } 2>&1 | awk '{
#    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
#           $1/3600, $1%3600/60, $1%60
#}' | tail -1
#echo ""


# filenames for scaffolder output
poa_file="./tmp/poa.fasta"
gr_file="./tmp/gr.fasta"
gm_file="./tmp/gm.fasta"
#poa_ext_file="./tmp/poa_ext.fasta"
poa_ext_file="$poa_dir/./extensions.fasta"
#gr_ext_file="./tmp/gr_ext.fasta"
gr_ext_file="$bwa_dir/./extensions.fasta"
#gm_ext_file="./tmp/gm_ext.fasta"
gm_ext_file="$graphmap_dir/./extensions.fasta"
poa_aln_file="./tmp/poa.sam"
gr_aln_file="./tmp/gr.sam"
gm_aln_file="./tmp/gm.sam"


# running scaffolder with global realign
echo "[EAGLER] extending contigs using global realignment method and bwa. Writing results to $bwa_dir"
echo "..."
{ time ./release/$name $2 $1 $bwa_dir; } 2>&1 | awk '{
    printf "[EAGLER] extension finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


# running scaffolder with poa
echo "[EAGLER] extending contigs using POA consensus method and bwa. Writing results to $poa_dir"
echo "..."
{ time ./release/$name -p $2 $1 $poa_dir; } 2>&1 | awk '{
    printf "[EAGLER] extension finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


# running scaffolder with graphmap
echo "[EAGLER] extending contigs using global realignment method and GraphMap. Writing results to $graphmap_dir"
echo "..."
{ time ./release/$name -g $2 $1 $graphmap_dir; } 2>&1 | awk '{
    printf "[EAGLER] extension finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "-----------------------------------------------------------------"
echo "Preparing results for analysis"
echo "-----------------------------------------------------------------"
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
bwa_mem="bwa mem -t $num_threads -x pacbio $3 $gr_ext_file"
{ time $bwa_mem > $gr_aln_file; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "[BWA] aligning POA extended contigs to reference genome"
echo "..."
bwa_mem="bwa mem -t $num_threads -x pacbio $3 $poa_ext_file"
{ time $bwa_mem > $poa_aln_file; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "[BWA] aligning GraphMap extended contigs to reference genome"
echo "..."
bwa_mem="bwa mem -t $num_threads -x pacbio $3 $gm_ext_file"
{ time $bwa_mem > $gm_aln_file; } 2>&1 | awk '{
    printf "[BWA] alignment finished in %d hours, %d minutes and %.3f seconds\n",
           $1/3600, $1%3600/60, $1%60
}' | tail -1
echo ""


echo "[ANALYSIS] analyzing .sam files"
python3 scripts/extension_analysis.py $3 $poa_aln_file $gr_aln_file $gm_aln_file
echo "[ANALYSIS] analysis finshed"

rm -rf "./tmp"
