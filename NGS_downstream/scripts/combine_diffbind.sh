#!/bin/bash
# Combines consensus and differential BEDs and prepares files for visualization in UCSC

input_dir=$1 # Contains BED files that were generated using /NGS_downstream/diffbind.R
output_dir=$2

cd $input_dir

for f in $(ls *_diff.bed | rev | cut -c 10- | rev | uniq)
do
  echo -e "Removes differential peaks from $f\_consensus.bed"
  bedtools intersect -v -a $f\_consensus.bed -b $f\_diff.bed > $f\_undiff.bed

  echo -e "Colors differential peaks for $f"
  awk 'BEGIN {OFS="\t"}; {if($9>0) $9="255,000,000"; else $9="000,000,255"; print $0;}' $f\_diff.bed > $f\_pre_rgb.bed
  awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4=".", $5="1000", $6=".", $7=$2, $8=$3, $9=$9}' $f\_pre_rgb.bed > $f\_rgb_diff.bed

  echo -e "Colors consensus peaks for $f"
  awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4=".", $5="1000", $6=".", $7=$2, $8=$3, $9="189,188,188"}' $f\_undiff.bed > $f\_rgb_undiff.bed

  echo -e "Combines differential and consensus peaks for $f"
  cat $f\_rgb_diff.bed $f\_rgb_undiff.bed > $f\_rgb_comb.bed
  grep 'chrX' $f\_rgb_comb.bed > $f\_rgb_xic.bed
  bedtools sort -i $f\_rgb_xic.bed > ${output_dir}$f\_rgb_sorted.bed
done
