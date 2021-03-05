#!/bin/bash
# Performs alignment and filtering of ChIP-seq data from Buecker et al. 2014 and Stadler et al. 2011

path=$1
fastq_dir=$2 # FASTQ files can be found at GSE56098 and GSE30203
work_dir=$3
ebwt=$4 # mm10 genome (prepared with build_bowtie2.sh)
echo "$ebwt"

cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd $fastq_dir

for f in $(ls *.fastq | rev | cut -c 7- | rev | uniq)
do
	cd ${data_dir}

	echo -e "Trimming $f with trim_galore"
	prun python3 trim_galore  --illumina ${fastq_dir}$f\.fastq >> $f\_trimmingStats.txt 2>&1

	echo -e "Mapping $f with bowtie2"
	bowtie2 --very-sensitive -x $ewbt -U $f\_trimmed.fq -S $f\.sam > $f\_mappingReport.txt
	rm $f\_trimmed.fq

	echo -e "Creating BAM files for $f and filtering for paired and mapped reads"
	samtools view -b -h -F 4 -q 20 $f\.sam > $f\.bam
	rm $f\.sam

	echo -e "Sorting BAM files for $f"
	samtools sort -m 1G $f\.bam -T $f\_sorted -o $f\_sorted.bam

	echo -e "Removing Blacklisted regions from $f"
	bedtools intersect -v -a $f\_sorted.bam -b $blacklistmm10 > $f\_sorted_blacklisted.bam

	echo -e "Removing duplicates from $f using PICARD\n"
	java -jar picard.jar MarkDuplicates INPUT=$f\_sorted_blacklisted.bam \
			OUTPUT=${bam_dir}$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT \
			REMOVE_DUPLICATES=TRUE
	samtools index ${bam_dir}$f\_dedup.bam
done
