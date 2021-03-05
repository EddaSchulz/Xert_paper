#Calculates and plots TPM in embryo rna-seq data from Zhang et al. 2018 for Fig. 4i and Supplementary Fig. 4f
library(tidyverse)
library(Rsubread)
library(Rgb)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] #Supply path to directory containing BAM files from Zhang et al. 2018 (generated via /NGS_alignment/RNAseq_align.sh)
gencode_xert <- args[2] #Supply path to GTF file containing GENCODE mm10 annotation with Xert added (At /NGS_downstream/files/GENCODE_vM25_plus_Xert.gtf)
output_dir <- args[3]

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

setwd(input_dir)
temp = list.files(pattern="*_dedup.bam$") #Files should be named after the following pattern: RNA_[Timepoint]_[Tissue]_[Replicate]_dedup.bam

feature_counts_embryo <- featureCounts(temp,
                                annot.ext = gencode_xert, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,  
                                GTF.featureType = "exon")
counts_embryo <- data.frame(feature_counts_embryo$counts, feature_counts_embryo$annotation) %>%
  select(-Chr, -Start, -End, -Strand)
names(counts_embryo) <- gsub(".dedup.bam", "", names(counts_embryo))
names(counts_embryo) <- gsub("RNA.", "", names(counts_embryo))

setwd(output_dir)
write_delim(counts_embryo, "Zhang_2018_counts.txt", delim = "\t")

TPM_embryo <- inner_join(counts_embryo, gene_names) %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length, -gene), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length, -GeneID)

write_delim(TPM_embryo, "Zhang_2018_TPM.txt", delim = "\t")
