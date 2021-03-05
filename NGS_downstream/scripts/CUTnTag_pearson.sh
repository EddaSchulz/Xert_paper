#!/bin/bash
#Pipeline to calculate Pearson correlation coefficient between CUT&Tag Samples in Gjaltema, Schwämmle et al. 2021

bam_dir=$1 #Contains BAM files generated using /NGS_alignment/CUTnTag_align.sh
merged_bam_dir=$2 #Contains merged BAM files generated using /NGS_downstream/merge_BAM.sh
chrM_bed=$3
output_dir=$4

cd $bam_dir

bamfiles=$(find -type f -regex '.*.bam' -printf '%f ')

echo -e "multiBamSummary for replicate CUT&Tag data"
prun python3 multiBamSummary bins -bs 1000 --bamfiles $bamfiles --smartLabels -bl $chrM_bed \
    -o ${output_dir}replicate_cnt.npz

echo -e "plotCorrelation for replicate CUT&Tag data"
prun python3 plotCorrelation --plotNumbers --outFileCorMatrix ${output_dir}replicate_cnt_pearson.txt \
    -in ${output_dir}replicate_cnt.npz -c pearson -p heatmap -o ${output_dir}deeptools_replicate_cnt_pearson.pdf


cd $merged_bam_dir

bamfiles=$(find -type f -regex '.*.bam' -printf '%f ')


echo -e "multiBamSummary for merged CUT&Tag data"
prun python3 multiBamSummary bins -bs 1000 --bamfiles $bamfiles --smartLabels -bl $chrM_bed \
    -o ${output_dir}merged_cnt.npz

echo -e "plotCorrelation for merged CUT&Tag data"
prun python3 plotCorrelation --plotNumbers --outFileCorMatrix ${output_dir}merged_cnt_pearson.txt \
    -in ${output_dir}merged_cnt.npz -c pearson -p heatmap -o ${output_dir}deeptools_merged_cnt_pearson.pdf
