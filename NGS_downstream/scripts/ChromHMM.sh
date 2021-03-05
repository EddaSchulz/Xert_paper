#!/bin/bash
# Script for annotating RE states using ChromHMM

path=$1
bam_dir=$2
output_dir=$3
binary_dir=$4


mm10=${path}software/ChromHMM/CHROMSIZES/mm10.txt
mark_table=${path}files/marktable.txt

echo -e "Creating binarized files"
java -mx4000M -jar ${path}software/ChromHMM/ChromHMM.jar BinarizeBam -paired -p 10 $mm10 $bam_dir $mark_table $binary_dir


echo -e "Learning model with 12 states"
java -mx4000M -jar ${path}software/ChromHMM/ChromHMM.jar LearnModel -p 10 $binary_dir $output_dir 12 mm10
