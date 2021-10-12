#!/bin/bash
# Master script to plot TPM in embryo rna-seq data from Zhang et al., 2018 and Deng et al., 2014 for Figures 3I and S3D
bam_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Calculates TPM from embryo RNA-seq data (Zhang et al., 2018/Deng et al., 2014) and plots Figure 3I and S3D"
   echo
   echo "Syntax: ./master_Embryo_Xert.sh [-b|d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "b     Provide directory containing processed BAM files from Zhang et al., 2018 RNA-seq. [mandatory]"
   echo "c     Provide directory containing processed BAM files from Deng et al., 2014 pseudobulk scRNA-seq. [mandatory]"
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -b)
            zhang_dir="$2"
            ;;
        -c)
            deng_dir="$2"
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

if [[ $zhang_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Zhang et al., 2018 BAM files with -b"
  exit 1
fi

if [[ $deng_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Deng et al., 2014 BAM files with -c"
  exit 1
fi

zhang_dir=$(realpath $zhang_dir)'/'
deng_dir=$(realpath $deng_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

gencode=${path}files/GENCODE_vM25_plus_Xert.gtf

# Counts reads and calculates TPM in different embryonic tissues
echo -e "Counting reads and calculating TPM"
Rscript ${path}scripts/TPM_Embryo_RNAseq.R $zhang_dir $deng_dir $gencode $work_dir

TPM_mat=${work_dir}Embryo_RNAseq_TPM.txt

# Plots Figures 3I and S3D
echo -e "Plotting expression data from embryo RNA-seq"
Rscript ${path}scripts/plot_Embryo_TPM.R $TPM_mat $work_dir
