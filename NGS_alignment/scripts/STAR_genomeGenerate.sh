#!/bin/bash
# Builds genome for STAR

path=$1
work_dir=$2
genome=$3 # hg38, mm10 and B6/Cast masked genome is found at /NGS_alignment/files/
gtf=$4 # hg38 and mm10 gene annotations

mkdir -p  ${work_dir}genome

STAR --runMode genomeGenerate --genomeDir ${work_dir}genome --genomeFastaFiles $genome --sjdbGTFfile $gtf
