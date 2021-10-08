#!/bin/bash
# Merges single cells according to their develpmental stage for Petropulos et al., 2016

path=$1
work_dir=$2
fastq_dir=$3 # Provide path to FASTQ files

cd $work_dir

mkdir -p pseudobulk_fastq
pseudo_fastq_dir=${work_dir}pseudobulk_fastq'/'

stage_info=${path}files/Petropoulos_sample_info.txt

cd ${fastq_dir}

while read old new; do
	echo -e "Copying $new"
  mv ${fastq_dir}$old\.fastq ${fastq_dir}$new
done < $stage_info

for n in E3 E4 E5 E6 E7;
do
  files=$(ls $n* | tr "\n" " ")
  echo -e "merging $n"
  cat $files > ${pseudo_fastq_dir}$n\_pseudobulk.fastq
done
