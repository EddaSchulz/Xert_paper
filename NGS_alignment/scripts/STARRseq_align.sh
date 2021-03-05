#!/bin/bash
# Performs alignment and filtering of STARR-seq data

path=$1
fastq_dir=$2 # FASTQ files can be found at GSEXXX
work_dir=$3
ebwt=$4 # B6/Cast masked genome (prepared with build_bowtie.sh)


cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd ${fastq_dir}

for f in $(ls *.fastq | rev | cut -c 9- | rev | uniq)
do
	cd ${data_dir}

	echo -e "Mapping $f with bowtie"
	bowtie -S -t -v 3 -m 1 -I 250 -X 2000 $ebwt -1 ${fastq_dir}$f\_1.fastq -2 ${fastq_dir}$f\_2.fastq $f\.sam >> $f\.mappingStats.txt 2>&1

	echo -e "Creating BAM files for $f and filtering for paired and mapped reads"
	samtools view -b -h -f 2 -q 10 $f\.sam > $f\.bam
	rm $f\.sam

	echo -e "Sorting BAM files for $f"
	samtools sort -m 1G $f\.bam -T $f\_sorted -o $f\_sorted.bam

	echo -e "Removing duplicates from $f using PICARD\n"
	java -jar picard.jar MarkDuplicates INPUT=$f\_sorted.bam \
			OUTPUT=${bam_dir}$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT \
	 		REMOVE_DUPLICATES=TRUE
	samtools index ${bam_dir}$f\_dedup.bam
done
