#!/bin/bash
# Master script to align scRNA-seq data from Petropoulos et al., 2016 as pseudobulk (E-MTAB-3929).
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Merges scRNAseq data from Petropoulos et al., 2016 per timepoint. Then builds hg38 genome and aligns the pseudobulk files."
   echo
   echo "Syntax: ./master_Petropoulos_pseudobulk_scRNAseq.sh [-d|f|p|h]"
   echo "options:"
   echo "d     Provide working directory. [optional]"
   echo "f     Provide directory containing FASTQ files (E-MTAB-3929). [mandatory]"
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
	echo -e "Please provide a path to FASTQ files from single cells from Petropoulos et al., 2016 (E-MTAB-3929) with -f"
  exit 1
fi

fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/hg38.fa
gtf=${path}files/GENCODE_vh38.gtf

# Merges single cells according to developmental timepoint
${path}scripts/Petropoulos_pseudobulk.sh $path $work_dir $fastq_dir

pseudobulk_fastq=${work_dir}pseudobulk_fastq'/'

# Build genome index for STAR from hg38 genome
${path}scripts/STAR_genomeGenerate.sh $path $work_dir $genome $gtf

ebwt=${work_dir}genome/

# Aligns FASTQ files and performs filtering steps
${path}scripts/Deng_Petropoulos_Soellner_RNAseq_align.sh $path $pseudobulk_fastq $work_dir $ebwt

# Merges replicate BAM files
${path}scripts/merge_stranded_BAM.sh $work_dir $bam_dir

final_bam_dir=${work_dir}final_bam'/'

# Generates coverage tracks
${path}scripts/generate_bigwig.sh $work_dir $final_bam_dir
