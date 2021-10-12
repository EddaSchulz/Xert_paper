#!/bin/bash
# Master script for data analysis and plotting based on MAGeCK
# Plots Figures 1D-F and SK-L
work_dir=$(pwd)'/'
path=''

help() {
   echo "Visualizes results from mageck mle analysis."
   echo
   echo "Syntax: ./master_mageck_analysis.sh [-d|p|h]"
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


# Plotting Figures 1D-E and S1K-L
echo -e "Plotting mageck mle vs genomic coordinates"
Rscript ${path}scripts/Xic_foldchange_plots.R ${path}files/ $work_dir

# Plotting Figure 1F
echo -e "Plotting betascore heatmap"
Rscript ${path}scripts/Betascore_heatmap.R ${path}files/ $work_dir
