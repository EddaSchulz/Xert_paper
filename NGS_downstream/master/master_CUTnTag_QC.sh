#!/bin/bash
# Master script to generate Supplementary Figures 2e-g
bam_dir=''
merged_bam_dir=''
work_dir=$(pwd)'/'
path=''
chip_dir=''
peak_dir=''

help() {
   echo "Generates and plots QC metrics for CUT&Tag data."
   echo
   echo "Syntax: ./master_CUTnTag_QC.sh [-b|d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "b     Provide directory containing replicate BAM files (from CUT&Tag). [mandatory]"
   echo "m     Provide directory containing merged BAM files (from CUT&Tag). [mandatory]"
   echo "n     Provide directory containing merged peak files (from CUT&Tag. [mandatory])"
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream. [mandatory]"
   echo "z     Provide directory containg BAM files from Zylicz et al. 2019. [mandatory]"
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
        -m)
            merged_bam_dir="$2"
            ;;
        -n)
            peak_dir="$2"
            ;;
        -z)
            chip_dir="$2"
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
	echo -e "Please provide the path to a directory containing replicate BAM files with -b"
  exit 1
fi

if [[ $merged_bam_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing merged BAM files with -m"
  exit 1
fi

if [[ $chip_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing ChIP-seq BAM files (Zylicz et al. 2019) with -z"
  exit 1
fi

if [[ $peak_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing merged peak files with -n"
  exit 1
fi


bam_dir=$(realpath $bam_dir)'/'
merged_bam_dir=$(realpath $merged_bam_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'
peak_dir=$(realpath $peak_dir)

chrM_bed=${path}files/chrM.bed

# Calculates Pearson correlation on replicates and merged BAM files for Supplementary Table 3
${path}scripts/CUTnTag_pearson.sh $bam_dir $merged_bam_dir $chrM_bed $work_dir

pearson_mat=${work_dir}merged_cnt_pearson.txt

# Plots correlation between samples as a heatmap
Rscript ${path}scripts/CUTnTag_pearson.R $pearson_mat $work_dir

# Prepares PCA analysis for CUT&Tag and native ChIP-seq
${path}scripts/nChIP_PCA.sh $chip_dir $merged_bam_dir $work_dir

pca_mat=${work_dir}nChIP_PCA.txt

#Plots PCA analaysis
Rscript ${path}scripts/nChIP_PCA.R $pca_mat $work_dir

#Performs binding enrichment using peakAnno
Rscript ${path}scripts/CUTnTag_peakAnno.R $peak_dir $work_dir
