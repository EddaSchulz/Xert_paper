#!/bin/bash
# Master script for read counting and DE gene enrichment for RNA-seq and TT-seq data
# Plots Fig. 4b and Supplementary Fig. 4b
tt_dir=''
rna_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Analyzes TT-seq and RNA-seq data to plot Fig. 4b and Sup. Fig. 4b."
   echo
   echo "Syntax: ./master_TTseq_downstream.sh [-d|h|p|r|t]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gjaltema_paper/NGS_downstream. [mandatory]"
   echo "r     Provide directory containing RNA-seq BAM files (/NGS_alignment/master/master_RNAseq.sh). [mandatory]"
   echo "t     Provide directory containing TT-seq BAM files (/NGS_alignment/master/master_TTseq.sh). [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -r)
            rna_dir="$2"
            ;;
        -t)
            tt_dir="$2"
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

if [[ $rna_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing RNA-seq BAM files with -r"
  exit 1
fi

if [[ $tt_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing TT-seq BAM files with -t"
  exit 1
fi

rna_dir=$(realpath $rna_dir)'/'
tt_dir=$(realpath $tt_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

gencode=${path}files/GENCODE_vM25_plus_Xert.gtf

# Counts reads, calculates TPM
echo -e "Counting reads and calculating TPM"
Rscript ${path}scripts/TPM_ttseq.R $tt_dir $rna_dir $gencode $work_dir


tt_counts=${work_dir}ttseq_counts.txt
tt_TPM=${work_dir}TPM_ttseq.txt
rna_counts=${work_dir}rnaseq_counts.txt
rna_TPM=${work_dir}TPM_rnaseq.txt

# Calculates DE genes using DEseq2
echo -e "Calculating DE genes between XXdXic and XO cell lines"
Rscript ${path}scripts/DEseq_ttseq.R $tt_counts $rna_counts $gencode $work_dir

DEseq_file=${work_dir}DEseq2_ttseq_total.txt

# Plots Xic gene expression in heatmap (Fig. 4b) and lineplot (Sup. Fig. 4b)
echo -e "Plots Xic lncRNA expression data"
Rscript ${path}scripts/plot_ttseq_TPM.R $tt_TPM $rna_TPM $DEseq_file $work_dir
