#!/bin/bash
# Script analyzes the Negative and High fractions of the CRISPR screen using MAGeCK test

count_table=$1 # Path to raw count table at /CRISPR_screen_analysis/files/raw_counts.txt
control_file=$2 # Path to list of non-targeting at /CRISPR_screen_analysis/files/nt_file.txt
output_dir=$3


prun python3 mageck test -k $count_table -t "R1_High,R2_High" -c "R1_Unsorted,R2_Unsorted" --norm-method control \
    --control-sgrna $control_file --output-prefix ${output_dir}Xic_screen
prun python3 mageck test -k $count_table -t "R1_Negative,R2_Negative" -c "R1_Unsorted,R2_Unsorted" --norm-method control \
    --control-sgrna $control_file --output-prefix ${output_dir}Xic_screen
