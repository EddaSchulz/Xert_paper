#This script calculates TPM for Figures 3I and S3D (Xert expression in embryonic tissues)
library(tidyverse)
library(Rsubread)
library(Rgb)

args <- commandArgs(trailingOnly = TRUE)
zhang_dir <- args[1] #Supply path to directory containing BAM files from Zhang et al., 2018 
deng_dir <- args[2] #Supply path to directory containing pseudobulk BAM files from Deng et al., 2014
gencode_xert <- args[3] #Supply path to GTF file containing GENCODE mm10 annotation with Xert added (At /NGS_downstream/files/GENCODE_vM25_plus_Xert.gtf)
output_dir <- args[4]

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

setwd(zhang_dir)
temp_zhang = list.files(pattern="*_dedup.bam$", full.names = TRUE) 

feature_counts_zhang <- featureCounts(temp_zhang,
                                annot.ext = gencode_xert, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,  
                                GTF.featureType = "exon")
counts_zhang <- data.frame(feature_counts_zhang$counts, feature_counts_zhang$annotation) %>%
  select(-Chr, -Start, -End, -Strand)

setwd(deng_dir)
temp_deng = list.files(pattern="*_dedup.bam$", full.names = TRUE) 

feature_counts_deng <- featureCounts(temp_deng,
                                      annot.ext = gencode_xert, isGTFAnnotationFile = TRUE, isPairedEnd = FALSE,  
                                      GTF.featureType = "exon")
counts_deng <- data.frame(feature_counts_deng$counts, feature_counts_deng$annotation) %>%
  select(-Chr, -Start, -End, -Strand)

counts <- inner_join(counts_deng, counts_zhang)
names(counts) <- gsub(x = names(counts), pattern = "\\.", replacement = "_")
names(counts) <- gsub(x = names(counts), pattern = "_5_", replacement = ".5_")
names(counts) <- gsub(x = names(counts), pattern = "_0_", replacement = ".0_")
names(counts) <- gsub(x = names(counts), pattern = "_cell", replacement = "cell")
names(counts) <- gsub(x = names(counts), pattern = "_dedup_bam", replacement = "")

setwd(output_dir)
TPM_embryo <- inner_join(counts, gene_names) %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length, -gene), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length, -GeneID) %>% 
  select(gene, everything())

write_delim(TPM_embryo, "Embryo_RNAseq_TPM.txt", delim = "\t")
