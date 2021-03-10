#!/bin/bash
# Master script for QC of the CRISPR screen
# Plots Supplementary Figures 1e and 1h-j
work_dir=$(pwd)'/'
path=''

help() {
   echo "Plots quality control for the CRISPR screen."
   echo
   echo "Syntax: ./master_screen_QC.sh [-d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/CRISPR_screen_analysis/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -d)
            work_dir="$2"
            ;;
        -p)
            path="$2"
            ;;
        -h)
            help
            exit 0
            ;;
        *)
            echo "Unknown or badly placed parameter '$1'." 1>&2
            exit 1
            ;;
    esac
}

while [[ "$#" -ge 1 ]]; do
    parse_args "$1" "$2"
    shift; shift
done

if [[ $path == '' ]]
then
	echo -e "Please provide the path to /Xert_paper/CRISPR_screen_analysis/ with -p"
  exit 1
fi



work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'


# Plotting Sup. Fig 1h
echo -e "Plotting correlation between replicates"
Rscript ${path}scripts/Screen_correlation.R ${path}files/ $work_dir

# Plotting Sup. Figs. 1e + 1i-j
echo -e "Plotting betascore heatmap"
Rscript ${path}scripts/guide_distribution.R ${path}files/ $work_dir
