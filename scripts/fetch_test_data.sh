#!/bin/bash

function print_delimiter {
    printf "%80s\n" | tr " " "="
}

if [ $# -eq 0 ]
then
    data_dir="./"
elif  [ $# -eq 1 ]
then
    data_dir=$1
    if [ ! -d ${data_dir} ]
    then
        echo "Destination directory does not exist"
        exit
    fi
else
    echo "usage: fetch_test_data.sh <dst_dir>"
    exit
fi

mkdir -p ${data_dir}

declare -A links

# all the tar files with nanopore reads
tars=(ERX708228.tar ERX708229.tar ERX708230.tar ERX708231.tar)

# links to data on EBI servers, usual DL speed 2-5 MB/s
links['ERX708228.tar']=ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar
links['ERX708229.tar']=ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar
links['ERX708230.tar']=ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar
links['ERX708231.tar']=ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar

print_delimiter
echo "Fetching missing files..."

# fetch all missing tar files
for (( i=0; i<${#tars[@]}; i++ ));
do
    tar_name=${tars[i]}
    tar_location=${data_dir}${tar_name}

    if [ -f ${tar_location} ]; then
        echo -e "\t${tar_name} found ($((${i}+1))/${#tars[@]})"
    else
        print_delimiter
        echo -e "\tfetching ${tar_name} ($((${i}+1))/${#tars[@]})"
        wget ${links[${tar_name}]} -O ${tar_location}
        print_delimiter
    fi
done

print_delimiter
