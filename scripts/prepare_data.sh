#!/bin/bash

function print_delimiter {
    printf "%80s\n" | tr " " "="
}

if [ $# -lt 1 ]
then
    echo "usage prepare_data.sh <dataset.tar>..."
    exit
fi

tars=("$@")

print_delimiter
echo "Preparing fast5 files..."

# prepare fast5 directories
for (( i=0; i<${#tars[@]}; i++ ));
do
    tar_id=$(basename -s .tar ${tars[${i}]})

    # create directory *.fast5/ directory for current dataset
    fast5_dir=${tar_id}\.fast5/
    mkdir -p ${fast5_dir}

    # extract the tar file
    echo -e "\t${tars[${i}]} extracting ($((${i}+1))/${#tars[@]})"
    tar -x -f ${tars[${i}]} -C ${fast5_dir}

    # group all fast5 files in the directory for the current dataset
    echo -e "\t${tars[${i}]} grouping fast5 files ($((${i}+1))/${#tars[@]})"
    find ${fast5_dir} -type f -name '*.fast5' -exec mv {} ${fast5_dir} \;

    # cleanup extracted folders
    find ${fast5_dir} -mindepth 1 -type d -delete
done

print_delimiter
echo -e "Extracting fasta reads..."

# extract reads with poretools
for (( i=0; i<${#tars[@]}; i++ ));
do
    tar_id=$(basename -s .tar ${tars[${i}]})
    fast5_dir=${tar_id}\.fast5/

    echo -e "\t${tar_id} extracting reads ($((${i}+1))/${#tars[@]})"
    poretools fasta --type 2D ${fast5_dir} >> raw.reads.unsorted
done

print_delimiter
