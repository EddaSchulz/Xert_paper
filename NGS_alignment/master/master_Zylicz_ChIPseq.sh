#!/bin/bash
# Master script to align ChIPseq data from Zylicz et al. 2019 (GSE116990)
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds masked mm10 genome, aligns and processes ChIPseq data from Zylicz et al. 2019."
   echo
   echo "Syntax: ./master_Zylicz_ChIPseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (GSE116990). [mandatory]"
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
	echo -e "Please provide a path to ChIPseq FASTQ files from Zylicz et al. (GSE116990) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/N_masked_B6_Cast.fa


# Build genome index for bowtie2 from N_masked mm10 genome
${path}scripts/build_bowtie2.sh $path $work_dir $genome

ebwt=${work_dir}genome/mm10

# Aligns FASTQ files and performs filtering steps
${path}scripts/CUTnTag_align.sh $path $fastq_dir $work_dir $ebwt

bam_dir=${work_dir}final_bam'/'

# Merges replicate BAM files
${path}scripts/merge_blacklisted_BAM.sh $work_dir $bam_dir
