#!/bin/bash
# Merges files ending on plus_dedup.bam and minus_dedup.bam

work_dir=$1
bam_dir=$2 # Contains BAM files

cd $work_dir

mkdir -p merged_bam
merged_bam_dir=${work_dir}merged_bam'/'

cd $bam_dir

for f in $(ls *_r[1-2]_dedup_plus.bam | rev | cut -c 19- | rev | uniq)
do
  echo -e "Merging BAM files for $f"
  samtools merge ${merged_bam_dir}$f\_merged_plus.bam $f\_r1_dedup_plus.bam $f\_r2_dedup_plus.bam
  samtools index ${merged_bam_dir}$f\_merged_plus.bam
done

for f in $(ls *_r[1-2]_dedup_minus.bam | rev | cut -c 20- | rev | uniq)
do
  echo -e "Merging BAM files for $f"
  samtools merge ${merged_bam_dir}$f\_merged_minus.bam $f\_r1_dedup_minus.bam $f\_r2_dedup_minus.bam
  samtools index ${merged_bam_dir}$f\_merged_minus.bam
done
