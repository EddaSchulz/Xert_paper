#!/bin/bash
# Merges single cells according to their develpmental stage for Deng et al., 2014

path=$1
work_dir=$2
fastq_dir=$3 # Provide path to FASTQ files

cd $work_dir

mkdir -p pseudobulk_fastq
pseudo_fastq_dir=${work_dir}pseudobulk_fastq'/'
stage_info=${path}files/Deng_sample_info.txt

cd ${fastq_dir}

while read srr stage; do
	echo -e "Changing names..."
  mv $srr\.fastq $srr\_$stage\.fastq
done < $stage_info

for f in $(ls *.fastq | rev | cut -c 7- | rev | grep -oP "_\K.*" | sort | uniq)
do
  echo -e "Merging $f"
  samples=$(ls *$f* | tr "\n" " ")
  cat $samples > ${pseudo_fastq_dir}$f\.fastq
done
