#!/bin/bash
# Master script to align RNA-seq data from Bauer et al., 2021 (GSE157448), Wang et al., 2019 (GSE116413) and Zhang et al. 2018 (GSE76505)
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds mm10 genome, aligns and processes published RNA-seq data from Bauer et al., 2021 (GSE157448), Wang et al., 2019 (GSE116413) or Zhang et al. 2018 (GSE76505)."
   echo
   echo "Syntax: ./master_Zhang_RNAseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (GSE157448/GSE116413/GSE76505). [mandatory]"
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_alignment/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -f)
            fastq_dir="$2"
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
	echo -e "Please provide a path to /Xert_paper/NGS_alignment/ with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide a path to RNA-seq FASTQ files from Bauer et al., 2021 (GSE157448), Wang et al., 2019 (GSE116413) or Zhang et al. 2018 (GSE76505) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/mm10.fa
gtf=${path}files/GENCODE_vM25_plus_Xert.gtf


# Build genome index for STAR from mm10 genome
${path}scripts/STAR_genomeGenerate.sh $path $work_dir $genome $gtf

ebwt=${work_dir}genome/

# Aligns FASTQ files and performs filtering steps
${path}scripts/RNAseq_align.sh $path $fastq_dir $work_dir $ebwt
