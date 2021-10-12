#!/bin/bash
# Master script for the LFC heatmap in Figures 2J
cnt_dir=''
atac_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Calculates LFC to plot heatmap in Figures 2J."
   echo
   echo "Syntax: ./master_RE_heatmap.sh [-a|b|c|d|p|t|h]"
   echo "options:"
   echo "a     Provide directory containing BAM files (Buecker et al., 2014). [mandatory]"
   echo "b     Provide directory containing BAM files (Wang et al., 2017). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            bueck_dir="$2"
            ;;
        -b)
            wang_dir="$2"
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
	echo -e "Please provide the path to /Xert_paper/NGS_downstream/ with -p"
  exit 1
fi

if [[ $bueck_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Buecker et al., 2014 BAM files with -a"
  exit 1
fi

if [[ $wang_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Wang et al., 2017 BAM files with -b"
  exit 1
fi

bueck_dir=$(realpath $bueck_dir)'/'
wang_dir=$(realpath $wang_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

# Plots heatmap
echo -e "Plots Figure 2J"
Rscript ${path}scripts/TF_heatmap.R $bueck_dir $wang_dir ${path}files/candidate_re.saf $work_dir
