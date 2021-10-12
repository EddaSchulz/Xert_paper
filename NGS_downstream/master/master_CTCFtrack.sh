#!/bin/bash
# Master script to identify CTCF binding sites for Figures 7A-B and S7A-B
work_dir=$(pwd)'/'
path=''

help() {
   echo "Prepares CTCF .tsv file from FIMO for use in UCSC (to visualize CBS)."
   echo
   echo "Syntax: ./master_CTCFtrack.sh [-d|p|h]"
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

# Uses .tsv file from FIMO to create colored BED track for UCSC
echo -e "Generating CBS track"
Rscript ${path}scripts/generate_CBS_bed.R ${path}files/ $work_dir
