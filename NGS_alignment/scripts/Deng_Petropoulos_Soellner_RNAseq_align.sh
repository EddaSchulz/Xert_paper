#!/bin/bash
# Performs alignment and filtering of single-end RNA-seq data

path=$1
fastq_dir=$2
work_dir=$3
ebwt=$4 # genome prepared with STAR_genomeGenerate.sh


cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd $fastq_dir

for f in $(ls *.fastq | rev | cut -c 7- | rev | uniq)
do
	cd ${data_dir}

  echo -e "Mapping $f with STAR"
  STAR --genomeDir $ebwt --readFilesIn ${fastq_dir}$f\.fastq \
				--outSAMtype BAM Unsorted --outFileNamePrefix $f\_ --outSAMattributes NH HI NM MD

  echo -e "Filtering $f for properly paired reads"
  samtools view -q 7 -F 4 -b $f\_Aligned.out.bam > $f\_filtered.bam

  echo -e "Sorting BAM files for $f"
  samtools sort -m 1G $f\_filtered.bam -T $f\_sorted -o $f\_sorted.bam

  echo -e "Removing duplicates from $f using PICARD\n"
  java -jar picard.jar MarkDuplicates INPUT=$f\_sorted_blacklisted.bam \
	 		OUTPUT=${bam_dir}$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT \
	 		REMOVE_DUPLICATES=TRUE
  samtools index ${bam_dir}$f\_dedup.bam
done
