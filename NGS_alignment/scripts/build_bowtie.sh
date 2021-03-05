#!/bin/bash
# Builds genome for bowtie

path=$1
work_dir=$2
genome=$3 # Standard mm10 is found at /NGS_alignment/files/

mkdir -p ${work_dir}genome

bowtie-build $genome ${work_dir}genome/mm10
