#!/bin/bash
# Calls peaks for ChIPseq data from Wang et al. 2107

work_dir=$1
bam_dir=$2 # Contains BAM files generated using /NGS_alignment/Wang_ChIPseq_align.sh

cd $work_dir

mkdir -p peaks
peaks_dir=${work_dir}peaks'/'

mkdir -p merged_peaks
merged_peaks_dir=${work_dir}merged_peaks'/'

cd $bam_dir

for f in $(ls *_dedup.bam | grep -v 'Input' | rev | cut -c 11- | rev | uniq)
do
  cond=$(echo "$f" | grep -oP '^.*?_\K(.*)')
  echo -e "Calling peaks for $f"
  macs2 callpeak -t ${bam_dir}$f\_dedup.bam -c ${bam_dir}Input_${cond}_dedup.bam -f BAMPE -g mm \
      -q 0.05 -n $f --outdir ${peaks_dir}
done
