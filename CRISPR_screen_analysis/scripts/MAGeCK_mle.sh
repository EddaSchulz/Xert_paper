#!/bin/bash
# Script analyzes the count table generated from MAGeCK count using MAGeCK mle

count_table=$1 # Path to raw count table at /CRISPR_screen_analysis/files/raw_counts.txt
design_matrix=$2 # Path to design matrix at /CRISPR_screen_analysis/files/design_matrix.txt
control_file=$3 # Path to list of non-targeting at /CRISPR_screen_analysis/files/nt_file.txt
output_dir=$4

prun python3 mageck mle -k $count_table -d $design_matrix --control-sgrna $control_file --norm-method control \
    --max-sgrnapergene-permutation 350 --output-prefix ${output_dir}Xic_screen
