#!/bin/bash
# Characterizes the Candidate REs based on which peaks were present in them
# Plots Figures S1B-D
work_dir=$(pwd)'/'
path=''

help() {
   echo "Characterizes candidate REs based on origin of peaks, gDNA number and length"
   echo
   echo "Syntax: ./master_characterize_REs.sh [-d|h|p]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/CRISPR_library_design/. [mandatory]"
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
	echo -e "Please provide the path to /Xert_paper/CRISPR_library_design/ with -p"
  exit 1
fi

work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

# Characterizes Candidate REs
echo -e "Characterizing candidate REs"
Rscript ${path}scripts/Char_Candidate_REs.R ${path}files $work_dir

# Plots Figures S1B-D
echo -e "Plotting candidate REs"
Rscript ${path}scripts/Plot_Candidate_REs.R ${work_dir}Char_candidate_RE.txt $work_dir
