#!/bin/bash
# Performs alignment and filtering of GRO-seq data from Wang et al., 2015 (GSE54471)

path=$1
fastq_dir=$2 # FASTQ files can be found at GSE54471
work_dir=$3
ebwt=$4 # hg38 genome (prepared with STAR_genomeGenerate.sh)


cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

mkdir -p strand_bam
strand_bam_dir=${work_dir}strand_bam'/'

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
  samtools sort -m 1G $f\_filtered.bam -T $f\_sorted -o ${bam_dir}$f\_sorted.bam

	echo -e "Splitting the BAM files according to the strand for $f"
	${path}scripts/split_allele.sh ${bam_dir}$f\_sorted.bam ${data_dir}
	mv $f\_sorted_minus.bam ${strand_bam_dir}$f\_sorted_minus.bam
	mv $f\_sorted_plus.bam ${strand_bam_dir}$f\_sorted_plus.bam
	samtools index ${strand_bam_dir}$f\_sorted_minus.bam
	samtools index ${strand_bam_dir}$f\_sorted_plus.bam
done
