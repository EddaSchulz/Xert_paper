#!/bin/bash
# Master script to design CRISPR screen library
# Returns BED files for Figure 1B
starr_dir=''
atac_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Combines ATAC-seq and STARR-seq data with FANTOM5 enhancers to create candidate REs withiin the Xic"
   echo
   echo "Syntax: ./master_library_design.sh [-d|h|p|s|a]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/CRISPR_library_design/. [mandatory]"
   echo "s     Provide directory containing STARR-seq BAM files (/NGS_alignment/master/master_STARRseq.sh). [mandatory]"
   echo "a     Provide directory containing ATAC-seq BAM files (/NGS_alignment/master/master_ATACseq.sh). [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            atac_dir="$2"
            ;;
        -s)
            starr_dir="$2"
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
	echo -e "Please provide the path to /Xert_paper/CRISPR_library_design/ with -p"
  exit 1
fi

if [[ $atac_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing ATAC-seq BAM files with -a"
  exit 1
fi

if [[ $starr_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing STARR-seq BAM files with -s"
  exit 1
fi

atac_dir=$(realpath $atac_dir)'/'
starr_dir=$(realpath $starr_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

mkdir -p ${work_dir}lib_design_files
lib_files=${work_dir}lib_design_files'/'

mkdir -p ${work_dir}lib_design_output
lib_out=${work_dir}lib_design_output'/'

# Returns peak file from STARR_seq
echo -e "Calling Peaks for STARR-seq"
#${path}scripts/STARRseq_peaks.sh $starr_dir $lib_files

# Returns peak file from ATAC_seq
echo -e "Calling Peaks for ATAC-seq"
${path}scripts/ATACseq_peaks.sh $atac_dir $lib_files

fantom5=${path}files/F5.mm10.enhancers.bed

# Limits Fantom5 enhancers to Xic
echo -e "Extracting Xic peaks from FANTOM5 data"
${path}scripts/FANTOM5_XIC.sh $fantom5 $lib_files

# Combines all three peak sets
echo -e "Combining all peaksets for Candidate RE set"
${path}scripts/Combine_Candidate_REs.sh ${lib_files}STARR_peaks_XIC.bed \
${lib_files}ATAC_peaks_XIC.bed ${lib_files}F5_enhancers_XIC.bed $lib_out


#At this point the Candidate_RE_merged.bed file was used to create the gRNA library using guidescan.com
#Afterwards, RE's longer than 2000 bps were split manually according to the ATAC-seq peaks
#Short RE's that were close together were combined (< 2000 bps combined size including the sequence in between)
