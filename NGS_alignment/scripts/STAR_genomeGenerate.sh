#!/bin/bash
# Builds genome for STAR

path=$1
work_dir=$2
genome=$3 # Standard mm10 and B6/Cast masked genome is found at /NGS_alignment/files/
gtf=$4

mkdir -p  ${work_dir}genome

STAR --runMode genomeGenerate --genomeDir ${work_dir}genome --genomeFastaFiles $genome --sjdbGTFfile $gtf
