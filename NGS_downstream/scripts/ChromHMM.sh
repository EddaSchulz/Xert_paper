#!/bin/bash
# Script for annotating RE states using ChromHMM

path=$1
atac_dir=$2
cnt_dir=$3
output_dir=$4
binary_dir=$5

mkdir -p ${output_dir}bam
bam_dir=${output_dir}bam'/'

cd $atac_dir
for f in $(ls *XX*)
do
  cp $f ${bam_dir}$f
done

cd $cnt_dir
for f in $(ls *XX*{H3K4me1,H3K4me3,H3K27ac,H3K27me3}*)
do
  cp $f ${bam_dir}$f
done

mm10=${path}software/ChromHMM/CHROMSIZES/mm10.txt
mark_table=${path}files/marktable.txt

echo -e "Creating binarized files"
java -mx4000M -jar ChromHMM.jar BinarizeBam -paired $mm10 $bam_dir $mark_table $binary_dir


echo -e "Learning model with 12 states"
java -mx4000M -jar ChromHMM.jar LearnModel $binary_dir $output_dir 12 mm10

rm -rf $bam_dir
