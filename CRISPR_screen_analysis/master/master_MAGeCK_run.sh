#!/bin/bash
# Master script to run MAGeCK count, mle and test
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Uses MAGeCK to align CRISPR screen data and perform read counting and statistical analysis."
   echo
   echo "Syntax: ./master_MaGeCK_run.sh [-f|d|p|h]"
   echo "options:"
   echo "f     Provide directory containing FASTQ files (CRISPR Screen). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/CRISPR_screen_analysis. [mandatory]"
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
	echo -e "Please provide the path to /Xert_paper/CRISPR_screen_analysis with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing CUT&Tag BAM files with -c"
  exit 1
fi


fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

mkdir -p ${work_dir}mageck_output
mageck_output=${work_dir}mageck_output'/'

# Aligns reads and performs counting using MAGeCK
echo -e "Counts reads using MAGeCK count"
${path}scripts/MAGeck_count.sh $fastq_dir ${path}files/library_file.txt ${path}files/nt_file.txt $work_dir

counts=${work_dir}Xic_screen.txt

# Performs statistical analysis using mageck mle
echo -e "Running MAGeCK mle"
${path}scripts/MAGeck_mle.sh $counts ${path}files/design_matrix.txt ${path}files/nt_file.txt $work_dir

# Performs statistical analysis using mageck test
echo -e "Running MAGeCK test"
${path}scripts/MAGeck_test.sh $counts ${path}files/nt_file.txt $work_dir
