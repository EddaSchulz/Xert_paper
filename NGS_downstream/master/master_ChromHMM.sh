#!/bin/bash
# Master script to create Fig. 2d and Supplementary Fig. 2h
atac_dir=''
cnt_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Uses ChromHMM to create Fig. 2d and Supplementary Fig. 2h."
   echo
   echo "Syntax: ./master_ChromHMM.sh [-b|d|p|h]"
   echo "options:"
   echo "a     Provide directory containing merged BAM files (ATAC-seq). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "c     Provide directory containing merged BAM files (CUT&Tag). [mandatory]"
   echo "h     Prints this help."
   echo "p     Provide path to /Xert_paper/NGS_downstream. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -a)
            atac_dir="$2"
            ;;
        -c)
            cnt_dir="$2"
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
	echo -e "Please provide the path to /Gjaltema_paper/NGS_downstream with -p"
  exit 1
fi

if [[ $atac_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing ATAC BAM files with -a"
  exit 1
fi

if [[ $cnt_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing CUT&Tag BAM files with -c"
  exit 1
fi


atac_dir=$(realpath $atac_dir)'/'
cnt_dir=$(realpath $cnt_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

mkdir -p ${work_dir}ChromHMM
chromm_dir=${work_dir}ChromHMM'/'

mkdir -p ${work_dir}Binarize_BAM
binary_dir=${work_dir}Binarize_BAM'/'

# Sections chromatin based on epigenetic marks for 12 states
${path}scripts/ChromHMM.sh $path $atac_dir $cnt_dir $chromm_dir $binary_dir


# After visual inspection of the resulting ${work_dir}ChromHMM/emissions_12.png file the states should be assigned to
# no RE, poised RE, weak RE and strong RE
echo -e "Enter 'strong RE' states as 'a b c' (space-delimited]"
read -a strongRE
printf -v strongRE_joined '%s,' "${strongRE[@]}"
strongRE_input=$(echo "${strongRE_joined%,}")

echo -e "Enter 'weak RE' states as 'a b c' (space-delimited]"
read -a weakRE
printf -v weakRE_joined '%s,' "${weakRE[@]}"
weakRE_input=$(echo "${weakRE_joined%,}")

echo -e "Enter 'poised RE' states as 'a b c' (space-delimited]"
read -a poisedRE
printf -v poisedRE_joined '%s,' "${poisedRE[@]}"
poisedRE_input=$(echo "${poisedRE_joined%,}")

echo -e "Enter 'no RE' states as 'a b c' (space-delimited]"
read -a noRE
printf -v noRE_joined '%s,' "${noRE[@]}"
noRE_input=$(echo "${noRE_joined%,}")

order=( "${noRE[@]}" "${poisedRE[@]}" "${weakRE[@]}" "${strongRE[@]}" )
printf -v order_joined '%s,' "${order[@]}"
order_input=$(echo "${order_joined%,}")

mkdir -p ${work_dir}Reduced_RE
red_dir=${work_dir}Reduced_RE'/'



# Plots a heatmap of the enrichment of the different data sets in the states
echo -e "Plotting emission state heatmap"
Rscript ${path}scripts/ChromHMM.R ${chromm_dir}emissions_12.txt $work_dir $order_input

cd $chromm_dir

# Reduces RE states in dense BED files
for f in $(ls *_dense.bed | rev | cut -c 11- | rev | uniq)
do
  echo -e "Reducing states for $f"
  ${path}scripts/reduce_states.py ${chromm_dir}$f\_dense.bed ${red_dir}$f\_reduced.bed $strongRE_input $weakRE_input $poisedRE_input $noRE_input
done
