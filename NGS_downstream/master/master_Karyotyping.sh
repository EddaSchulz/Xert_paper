#!/bin/bash
# Master script to plot karyotyping heatmaps for Figure S4F
work_dir=$(pwd)'/'
path=''

help() {
   echo "Prepares heatmap for karyotyping data of all novel cell lines."
   echo
   echo "Syntax: ./master_karyotypiing.sh [-d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream/. [mandatory]"
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
	echo -e "Please provide the path to /Xert_paper/NGS_downstream/ with -p"
  exit 1
fi

work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

# Uses karyotyping plots to generate a chromosome-seperated heatmap
echo -e "Generating karyotyping heatmap"
Rscript ${path}scripts/karyotyping_plot.R ${path}files/ $work_dir
