#!/bin/bash
# Master script to align GRO-seq data from Wang et al., 2015 (GSE54471)
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds hg38 genome, aligns and processes RNA-seq data."
   echo
   echo "Syntax: ./master_Wang_GROseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (GSE54471). [mandatory]"
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
	echo -e "Please provide a path to Wang et al., 2015 GRO-seq FASTQ files (GSE54471) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/hg38.fa
gtf=${path}files/GENCODE_vH38.gtf

# Build genome index for STAR from hg38 genome
${path}scripts/STAR_genomeGenerate.sh $path $work_dir $genome $gtf

ebwt=${work_dir}genome/

# Aligns FASTQ files and performs filtering steps
${path}scripts/Wang_GROseq_align.sh $path $fastq_dir $work_dir $ebwt

strand_bam_dir=${work_dir}strand_bam'/'

# Generates coverage tracks
${path}scripts/generate_bigwig.sh $work_dir $strand_bam_dir
