#!/bin/bash
# Master script to plot TPM in embryo rna-seq data from Zhang et al. 2018 for Fig. 4i and Supplementary Fig. 4f
bam_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Calculates TPM from embry RNA-seq data (Zhang et al. 2018) and plots Fig. 4i and Supplementary Fig. 4f."
   echo
   echo "Syntax: ./master_Zhang_RNAseq_downstream.sh [-b|d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "b     Provide directory containing processed BAM files (from NGS_alignment/master_Zhang_RNAseq.sh). [mandatory]"
   echo "h     Prints this help."
   echo "p     Provide path to /Gjaltema_paper/NGS_downstream. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -b)
            bam_dir="$2"
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
	echo -e "Please provide the path to /Gjaltema_paper/NGS_downstream with -p"
  exit 1
fi

if [[ $bam_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing BAM files with -b"
  exit 1
fi

bam_dir=$(realpath $bam_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

gencode=${path}files/GENCODE_vM25_plus_Xert.gtf

# Counts reads and calculates TPM in different embryonic tissues
echo -e "Counting reads and calculating TPM"
Rscript ${path}scripts/TPM_Zhang_2018.R $bam_dir $gencode $work_dir

TPM_mat=${work_dir}Zhang_2018_TPM.txt

# Plots Fig. 4i and Sup. Fig. 4f
echo -e "Plotting expression data from Zhang et al. 2018"
Rscript ${path}scripts/plot_Zhang_2018.R $TPM_mat $work_dir
