#!/bin/bash
# Master script for differential binding analysis using Diffbind and diffReps
cnt_dir=''
atac_dir=''
cnt_peaks=''
atac_peaks=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Performs differential binding analysis for CUT&Tag and ATAC-seq."
   echo
   echo "Syntax: ./master_diffEnrichment.sh [-a|b|c|d|p|t|h]"
   echo "options:"
   echo "a     Provide directory containing merged peak files (ATAC-seq). [mandatory]"
   echo "b     Provide directory containing BAM files (CUT&Tag. [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "c     Provide directory containing merged peak files (CUT&Tag). [mandatory]"
   echo "h     Prints this help."
   echo "p     Provide path to /Gjaltema_paper/NGS_downstream. [mandatory]"
   echo "t     Provide directory containing BAM files (ATAC-seq). [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            atac_peaks="$2"
            ;;
        -b)
            cnt_dir="$2"
            ;;
        -c)
            cnt_peaks="$2"
            ;;
        -d)
            work_dir="$2"
            ;;
        -p)
            path="$2"
            ;;
        -t)
            atac_dir="$2"
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
	echo -e "Please provide the path to /Xert_paper/NGS_downstream with -p"
  exit 1
fi

if [[ $cnt_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing CUT&Tag BAM files with -b"
  exit 1
fi

if [[ $atac_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing ATAC-seq BAM files with -t"
  exit 1
fi


if [[ $cnt_peaks == '' ]]
then
	echo -e "Please provide the path to a directory containing CUT&Tag peak files with -c"
  exit 1
fi

if [[ $atac_peaks == '' ]]
then
	echo -e "Please provide the path to a directory containing ATAC-seq peak files with -a"
  exit 1
fi

atac_peaks=$(realpath $atac_peaks)'/'
cnt_peaks=$(realpath $cnt_peaks)'/'
atac_dir=$(realpath $atac_dir)'/'
cnt_dir=$(realpath $cnt_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

mkdir -p ${work_dir}noX_peaks
noX_dir=${work_dir}noX_peaks'/'

# Prepares BED files without X chromosomes for Diffbind analysis
echo -e "Prepares BED files for DiffBind"
${path}scripts/diffbind.sh $atac_peaks $cnt_peaks $noX_dir

mkdir -p ${work_dir}raw_DiffBind
rawDiffbind_dir=${work_dir}raw_DiffBind'/'

# Performs Diffbind analysis in R
echo -e "Performs Diffbind analysis"
Rscript ${path}scripts/diffbind.R $noX_dir $atac_dir $cnt_dir $rawDiffbind_dir

mkdir -p ${work_dir}final_DiffBind
DiffBind_dir=${work_dir}final_DiffBind'/'

# Prepares BED files for visualization with UCSC
echo -e "Prepares BED tracks for UCSC"
${path}scripts/combine_diffbind.sh $rawDiffbind_dir $DiffBind_dir

chrom_sizes=${path}files/mm10_chrom_sizes.txt

mkdir -p ${work_dir}raw_diffreps
raw_diffreps=${work_dir}raw_diffreps'/'

# Runs diffreps to call differential regions between H3K9me3 samples
echo -e "Runs diffReps for H3K9me3"
${path}scripts/diffreps.sh $cnt_dir $raw_diffreps $chrom_sizes

# Runs diffreps to call consensus regions between H3K9me3 samples
echo -e "Calls consensus peaks between H3K9me3 conditions"
${path}scripts/diffreps_consensus.sh $raw_diffreps $chrom_sizes

mkdir -p ${work_dir}final_diffreps
diffreps_dir=${work_dir}final_diffreps'/'

# Combines consensus and differential peaks from diffreps
echo -e "Prepares diffReps BED tracks for UCSC"
${path}scripts/combine_diffreps.sh $raw_diffreps $diffreps_dir
