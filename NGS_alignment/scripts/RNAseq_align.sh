#!/bin/bash
# Performs alignment and filtering of RNA-seq data

path=$1
fastq_dir=$2 # FASTQ files can be found at GSE167356
work_dir=$3
ebwt=$4 # B6/Cast masked genome (prepared with STAR_genomeGenerate.sh)


cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd $fastq_dir

for f in $(ls *.fastq | rev | cut -c 9- | rev | uniq)
do
	cd ${data_dir}

  echo -e "Mapping $f with STAR"
  STAR --genomeDir $ebwt --readFilesIn ${fastq_dir}$f\_1.fastq ${fastq_dir}$f\_2.fastq \
				--outSAMtype BAM Unsorted --outFileNamePrefix $f\_ --outSAMattributes NH HI NM MD

  echo -e "Filtering $f for properly paired reads"
  samtools view -q 7 -f 3 -b $f\_Aligned.out.bam > $f\_filtered.bam

  echo -e "Sorting BAM files for $f"
  samtools sort -m 1G $f\_filtered.bam -T $f\_sorted -o $f\_sorted.bam

  echo -e "Removing blacklisted regions for $f"
  bedtools intersect -v -a $f\_sorted.bam -b ${path}files/mm10.bl.bed > $f\_sorted_blacklisted.bam

  echo -e "Removing duplicates from $f using PICARD\n"
  java -jar picard.jar MarkDuplicates INPUT=$f\_sorted_blacklisted.bam \
	 		OUTPUT=${bam_dir}$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT \
	 		REMOVE_DUPLICATES=TRUE
  samtools index ${bam_dir}$f\_dedup.bam
done
