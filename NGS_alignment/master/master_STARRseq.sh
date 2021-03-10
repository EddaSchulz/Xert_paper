#!/bin/bash
# Master script to align STARR-seq data, call peaks and create coverage tracks
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds mm10 genome, aligns and processes STARR-seq data."
   echo
   echo "Syntax: ./master_STARRseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (GSE167355). [mandatory]"
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
	echo -e "Please provide a path to STARR-seq FASTQ files (GSEXXX) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/mm10.fa


# Build genome index for bowtie from N_masked mm10 genome
${path}scripts/build_bowtie.sh $path $work_dir $genome

ebwt=${work_dir}genome/mm10

# Aligns FASTQ files and performs filtering steps
${path}scripts/STARRseq_align.sh $path $fastq_dir $work_dir $ebwt
