#!/bin/bash
# Master script to analyze CRISPR screen with a bootstrapping approach and plot Figures 1G and S1M
work_dir=$(pwd)'/'
path=''

help() {
   echo "Uses Bootstrapping approach to analyze CRISPRi screen data."
   echo
   echo "Syntax: ./master_bootstrap_analysis.sh [-d|p|h]"
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


# Running the script
echo -e "Performing statistical analysis using Bootstrapping approach"
Rscript ${path}scripts/bootstrapping_analysis.R ${path}files/raw_counts.txt $work_dir
