#!/bin/bash
# Script generates a candidate RE set within the XIC from the STARR-seq, ATAC-seq and FANTOM5 data

starr_peaks=$1
atac_peaks=$2
F5_enhancers=$3
output_dir=$4

cd $output_dir

cat $starr_peaks $atac_peaks $F5_enhancers > Candidate_RE_raw.bed
bedtools sort -i Candidate_RE_raw.bed > Candidate_RE_sorted.bed
bedtools merge -i Candidate_RE_sorted.bed > Candidate_RE_merged.bed

rm Candidate_RE_raw.bed
rm Candidate_RE_sorted.bed
