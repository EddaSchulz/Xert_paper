#!/bin/bash
# Master script for z-score and LFC heatmaps for Fig. 3c-d
cnt_dir=''
atac_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Calculates zscore and LFC to plot heatmaps for Fig. 3c-d."
   echo
   echo "Syntax: ./master_RE_heatmap.sh [-a|b|c|d|p|t|h]"
   echo "options:"
   echo "a     Provide directory containing BAM files (ATAC-seq). [mandatory]"
   echo "b     Provide directory containing BAM files (CUT&Tag. [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            atac_dir="$2"
            ;;
        -b)
            cnt_dir="$2"
            ;;
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
	echo -e "Please provide the path to /Xert_paper/NGS_downstream with -p"
  exit 1
fi

if [[ $cnt_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing CUT&Tag BAM files with -b"
  exit 1
fi

if [[ $atac_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing ATAC-seq BAM files with -a"
  exit 1
fi

atac_dir=$(realpath $atac_dir)'/'
cnt_dir=$(realpath $cnt_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

# Calculates Z-score and LFC and plots heatmaps
echo -e "Plots Fig. 3c-d"
Rscript ${path}scripts/RE_cnt_heatmaps.R $atac_dir $cnt_dir ${path}files/candidate_re.saf $work_dir
