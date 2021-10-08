#!/bin/bash
# Builds genome for bowtie2

path=$1
work_dir=$2
genome=$3

mkdir -p  ${work_dir}genome

bowtie2-build $genome ${work_dir}genome/mm10
