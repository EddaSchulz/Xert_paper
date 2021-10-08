#!/bin/bash
# Builds hg38 genome for bowtie2

path=$1
work_dir=$2
genome=$3 # hg38 is found at /NGS_alignment/files/

mkdir -p  ${work_dir}genome

bowtie2-build $genome ${work_dir}genome/hg38
