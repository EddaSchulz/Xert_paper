#!/bin/bash
# Script generates a composite peak set of the ATAC-seq data within the XIC

input_dir=$1 # Contains BAM files generated using /NGS_alignment/master/master_ATACseq.sh
output_dir=$2


cd $input_dir

echo -e "Merging treatment files"
bam_files=$(find -type f -regex ".*_dedup.bam" -printf "%f ")
samtools merge ${output_dir}ATAC_all_merged.bam $bam_files

cd $output_dir

echo -e "Calling peaks for the merged ATAC-seq data"
macs2 callpeak -t ATAC_all_merged.bam -f BAMPE -g mm -q 0.1 -n ATAC_all_merged --outdir ${output_dir}

awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' ATAC_all_merged_peaks.narrowPeak > ATAC_all_merged.bed
echo -e "chrX\t103198658\t104058961" > XIC_coords.bed
bedtools intersect -a XIC_coords.bed -b ATAC_all_merged.bed > ATAC_peaks_XIC.bed
