#!/bin/bash
# Master script to plot the results of the Cistrome-DB toolkit analysis for Sup. Fig. 3d
work_dir=$(pwd)'/'
path=''

help() {
   echo "Plots results from Cistrome-DB toolkit analysis on Xert REs."
   echo
   echo "Syntax: ./master_cistrome_search.sh [-d|p|h]"
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

# Uses TSV output from cistrome DB-toolkit to plot the enrichment of the top 15 TF's within Xert
echo -e "Plots Cistrome DB toolkit analysis"
Rscript ${path}scripts/cistrome_plot.R ${path}files/ $work_dir
