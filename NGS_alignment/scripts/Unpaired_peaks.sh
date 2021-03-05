#!/bin/bash
# Calls peaks for unpaired ChIP-seq data

work_dir=$1
bam_dir=$2 # Contains BAM files generated using /NGS_alignment/Unpaired_ChIPseq_align.sh

cd $work_dir

mkdir -p peaks
peaks_dir=${work_dir}peaks'/'

mkdir -p merged_peaks
merged_peaks_dir=${work_dir}merged_peaks'/'

cd $bam_dir

for f in $(ls *_r[1-2]_dedup.bam | grep -v 'Input' | rev | cut -c 14- | rev | uniq)
do
  cond=$(echo "$f" | grep -oP '^.*?_\K(.*)')

  echo -e "Calling peaks for $f"
  macs2 callpeak -t ${bam_dir}*_r1_dedup.bam -c ${bam_dir}Input_${cond}_r1_dedup.bam -g mm -q 0.05 -n $f\_r1 --outdir ${peaks_dir}
  macs2 callpeak -t ${bam_dir}*_r2_dedup.bam -c ${bam_dir}Input_${cond}_r2_dedup.bam -g mm -q 0.05 -n $f\_r2 --outdir ${peaks_dir}

  echo -e "Merge peak files for $f"
  bedtools intersect -a ${peaks_dir}$f\_r1_peaks.narrowPeak -b ${peaks_dir}$f\_r2_peaks.narrowPeak \
      > ${merged_peaks_dir}$f\_merged.narrowPeak
  awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' ${merged_peaks_dir}$f\_merged.narrowPeak \
      > ${merged_peaks_dir}$f\_merged.bed
  rm ${merged_peaks_dir}$f\_merged.narrowPeak
done
