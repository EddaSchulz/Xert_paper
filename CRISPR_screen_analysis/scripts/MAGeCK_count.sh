#!/bin/bash
# Script generates a count table from FASTQ files of the CRISPRi screen in Gjaltema, Schwämmle et al., 2021

fastq_dir=$1 # FASTQ files are found at GSE167352
library_file=$2 # Path to sgRNA library file (Can be found at /CRISPR_screen_analysis/files/library_file.txt)
control_file=$3 # Path to list of non-targeting controls (Can be found at /CRISPR_screen_analysis/files/nt_file.txt)
output_dir=$4


fastq_files=$(find -type f -regex "*.fastq" -printf "%f ")

prun python3 mageck count -l $library_file --fastq $fastq_files --norm-method control --control-sgrna $control_file \
    --output-prefix ${output_dir}Xic_screen
