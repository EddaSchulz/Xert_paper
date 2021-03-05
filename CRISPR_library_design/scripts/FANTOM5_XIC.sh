#!/bin/bash
# Script filters Fantom5 mm10 enhancers (Lizio et al. 2015, 2019) for the XIC

fantom5_enhancers=$1 # Path to Fantom5 mm10 enhancers (Lizio et al. 2015, 2019). File can be found at CRISPR_library_design/files/.
output_dir=$2

cd $output_dir

echo -e "Filtering FANTOM5 enhancer for the XIC"
echo -e "chrX\t103198658\t104058961" > XIC_coords.bed
bedtools intersect -a XIC_coords.bed -b $fantom5_enhancers > F5_enhancers_XIC.bed
