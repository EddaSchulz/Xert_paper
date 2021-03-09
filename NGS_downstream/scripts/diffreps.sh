#!/bin/bash
# Script to create differential enrichment tracks for H3K9me3 using diffreps

input_dir=$1 # Contains X_blacklisted_sorted.bam files for H3K9me3 that were generated using /NGS_alignment/CUTnTag_align.sh
output_dir=$2
chrom_sizes=$3 # Can be found at /NGS_downstream/files/mm10_chrom_sizes.txt

cd $input_dir

for f in $(ls *H3K9me3*_sorted_blacklisted.bam | rev | cut -c 24- | rev | uniq)
do
  echo -e "Creates BED file for $f"
  bedtools bamtobed -i $f\_sorted_blacklisted.bam > ${output_dir}$f\_pre.bed
done

cd $output_dir

echo -e "XO_H3K9me3_d0\tXXdXic_H3K9me3_d0\tclone_d0\nXO_H3K9me3_d2\tXXdXic_H3K9me3_d2\tclone_d2\nXO_H3K9me3_d4\tXXdXic_H3K9me3_d4\tclone_d4\nXXdXic_H3K9me3_d0\tXXdXic_H3K9me3_d2\tdiff_d2\nXXdXic_H3K9me3_d0\tXXdXic_H3K9me3_d4\tdiff_d4" > comps.txt

while read line
do
  sampleA=$(echo "$line" | cut -f1)
  sampleB=$(echo "$line" | cut -f2)
  comp=$(echo "$line" | cut -f3)
  echo -e "Conducts differential peak binding analysis using diffreps for $comp"
  diffReps.pl -tr CUTnTag_TX1072_$sampleA\_r1_pre.bed CUTnTag_TX1072_$sampleA\_r2_pre.bed \
      -co CUTnTag_TX1072_$sampleB\_r1_pre.bed CUTnTag_TX1072_$sampleB\_r2_pre.bed \
      -re $comp\_H3K9me3_diff.bed -ch $chrom_sizes  --window 5000 --step 1000 -me nb
done < comps.txt

rm comps.txt
