#!/bin/bash
# Master script to calculate TPM values for adult/differentiated RNA-seq samples (Soellner et al., 2017, Bauer et al., 2021, Wang et al., 2019)
cnt_dir=''
atac_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Calculates LFC to plot heatmap in Figures 2J."
   echo
   echo "Syntax: ./master_RE_heatmap.sh [-a|b|c|d|p|t|h]"
   echo "options:"
   echo "a     Provide directory containing BAM files from adult mouse tissues (Soellner et al., 2017). [mandatory]"
   echo "b     Provide directory containing BAM files from mouse NPCs (Bauer et al., 2021). [mandatory]"
   echo "c     Provide directory containing BAM files from mouse MEFs (Wang et al., 2019). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            soellner_dir="$2"
            ;;
        -b)
            bauer_dir="$2"
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

if [[ $soellner_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Soellner al., 2017 BAM files with -a"
  exit 1
fi

if [[ $bauer_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Bauer et al., 2021 BAM files with -b"
  exit 1
fi

if [[ $wang_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Wang et al., 2019 BAM files with -c"
  exit 1
fi

soellner_dir=$(realpath $soellner_dir)'/'
bauer_dir=$(realpath $bauer_dir)'/'
wang_dir=$(realpath $wang_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

gencode=${path}files/GENCODE_vM25_plus_Xert.gtf

# Calculates TPM and plots values for Xert
echo -e "Calculates TPM and plots figure S3E"
Rscript ${path}scripts/TPM_adult_Xert.R $soellner_dir $bauer_dir $wang_dir $gencode $work_dir
