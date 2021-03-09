#!/bin/bash
# Master script to align pA-RNAseq data and create coverage tracks
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds masked mm10 genome, aligns and processes pA-RNAseq data and creates coverage tracks."
   echo
   echo "Syntax: ./master_pA-RNAseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (GSE167354). [mandatory]"
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
	echo -e "Please provide a path to /Gjaltema_paper/NGS_alignment with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide a path to pA-RNAseq FASTQ files (GSE167354) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/N_masked_B6_Cast.fa
gtf=${path}files/GENCODE_vM25_plus_Xert.gtf


# Build genome index for STAR from N_masked mm10 genome
${path}scripts/STAR_genomeGenerate.sh $path $work_dir $genome $gtf

ebwt=${work_dir}genome/

# Aligns FASTQ files and performs filtering steps
${path}scripts/TTseq_align.sh $path $fastq_dir $work_dir $ebwt

bam_dir=${work_dir}strand_bam'/'

# Generates coverage tracks
${path}scripts/generate_bigwig.sh $work_dir $bam_dir
