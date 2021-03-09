#!/bin/bash
# Script to create consensus peak files for H3K9me3 using diffreps

work_dir=$1 # Contains BED files for H3K9me3 that were generated using /NGS_downstream/diffreps.sh
chrom_sizes=$2 # Can be found at /NGS_downstream/files/mm10_chrom_sizes.txt

cd $work_dir

for f in $(ls *_pre.bed | rev | cut -c 9- | rev | uniq)
do
  echo -e "Create BED with empty chrX for $f"
  sed '/chrX/d' $f\_pre.bed > $f\_noX.bed
done

for f in $(ls *_noX.bed | rev | cut -c 12- | rev | uniq)
do
  echo -e "Simulating peakcalling against empty chrX for $f"
  diffReps.pl -tr $f\_r1_pre.bed  $f\_r2_pre.bed -co ${output_dir}$f\_r1_noX.bed ${output_dir}$f\_r2_noX.bed \
    -re ${output_dir}$f\_diffreps.bed -ch $chrom_sizes --window 5000 --step 1000 -me nb
  sed -i '/chrX/!d' $f\_diffreps.bed
done

echo -e "XO_H3K9me3_d0\tXXdXic_H3K9me3_d0\tclone_d0\nXO_H3K9me3_d2\tXXdXic_H3K9me3_d2\tclone_d2\nXO_H3K9me3_d4\tXXdXic_H3K9me3_d4\tclone_d4\nXXdXic_H3K9me3_d0\tXXdXic_H3K9me3_d2\tdiff_d2\nXXdXic_H3K9me3_d0\tXXdXic_H3K9me3_d4\tdiff_d4" > comps.txt

while read line
do
  sampleA=$(echo "$line" | cut -f1)
  sampleB=$(echo "$line" | cut -f2)
  comp=$(echo "$line" | cut -f3)
  echo -e "Creates consensus peak file for $comp"
  bedtools intersect -a CUTnTag_TX1072_$sampleA\_diffreps.bed -b CUTnTag_TX1072_$sampleB\_diffreps.bed  \
      > $comp\_H3K9me3_consensus_unfiltered.bed
  echo -e "Removes differential peaks from consensus files for $comp"
  bedtools intersect -v -a $comp\_H3K9me3_consensus_unfiltered.bed -b $comp\_H3K9me3_diff.bed  > $comp\_H3K9me3_consensus.bed
done < comps.txt

rm comps.txt
