#!/bin/bash
# Prepares peak files for differential binding analysis with Diffbind

atac_dir=$1 # Contains merged peak files for ATAC-seq
cnt_dir=$2 # Contains merged peak files for CUT&Tag
output_dir=$3

echo -e "chrX\t1\t103182700\nchrX\t103955532\t171031299" > ${output_dir}chrX_minus_xic.bed

cd $atac_dir

for f in $(ls *_merged.bed | rev | cut -c 12- | rev | uniq)
do
  echo -e "Removes peaks on chrX outside Xic deletion for $f"
  bedtools intersect -v -a $f\_merged.bed -b ${output_dir}chrX_minus_xic.bed > ${output_dir}$f\_r1_noX.bed
  bedtools intersect -v -a $f\_merged.bed -b ${output_dir}chrX_minus_xic.bed > ${output_dir}$f\_r2_noX.bed
done

cd $cnt_dir

for f in $(ls *{H3K4me1,H3K4me3,H3K27ac}*_merged.bed | rev | cut -c 12- | rev | uniq)
do
  echo -e "Removes peaks on chrX outside Xic deletion for $f"
  bedtools intersect -v -a $f\_merged.bed -b ${output_dir}chrX_minus_xic.bed > ${output_dir}$f\_r1_noX.bed
  bedtools intersect -v -a $f\_merged.bed -b ${output_dir}chrX_minus_xic.bed > ${output_dir}$f\_r2_noX.bed
done


rm ${output_dir}chrX_minus_xic.bed
