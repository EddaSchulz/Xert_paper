#!/bin/bash
# Script generates a composite peak set of the STARR-seq data within the XIC

bam_dir=$1 # Contains BAM files generated using /NGS_alignment/master/master_STARRseq.sh
output_dir=$2


cd $bam_dir

echo -e "Merging treatment files"
bam_files=$(find -type f -regex ".*r[1-3]_dedup.bam" -printf "%f ")
samtools merge ${output_dir}STARR_1.8_all_merged.bam $bam_files


cd $output_dir

echo -e "Calling peaks for the merged STARR-seq data"
macs2 callpeak -t STARR_1.8_all_merged.bam -c ${bam_dir}STARR_1.8_Input_dedup.bam \
    -f BAMPE -g mm -q 0.1 -n STARR_all_merged --outdir ${output_dir}

awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' STARR_all_merged_peaks.narrowPeak > STARR_all_merged.bed
echo -e "chrX\t103198658\t104058961" > XIC_coords.bed
bedtools intersect -a XIC_coords.bed -b STARR_all_merged.bed > STARR_peaks_XIC.bed
