#!/bin/bash
# Prepares samples for PCA analysis comparing nChIP-seq and CUT&Tag data
nChip_dir=$1 # Contains merged BAM files from nChIP in TX1072 (Zylicz et al. 2019)
CnT_dir=$2 # Contains merged BAM files from CUT&Tag
output_dir=$3

cd $nChip_dir
nChip_files=$(find -type f -regex '.*.bam' | xargs realpath | xargs)

cd $CnT_dir
CnT_files=$(find -type f -regex '.*Xic_\(H3K4me1\|H3K4me3\|H3K27ac\|H3K27me3\|H2AK119ub\)_d0.*.bam' | xargs realpath | xargs)

echo -e "Running multiBamSummary"
prun python3 multiBamSummary bins -bs 1000 --bamfiles $nChip_files $CnT_files --smartLabels \
 -o ${output_dir}nChIP_PCA.npz --outRawCounts ${output_dir}nChIP_PCA.txt
