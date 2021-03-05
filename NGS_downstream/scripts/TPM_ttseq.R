#This script calculates TPM for Supplementary Table 4
library(Rsubread)
library(EnvStats)
library(Rgb)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
TT_dir <- args[1] #Supply path to directory containing TT-seq BAM files
RNA_dir <- args[2] #Supply path to directory containing RNA-seq BAM files 
gencode_xert <- args[3] #Supply path to GTF file containing GENCODE mm10 annotation with Xert added 
output_dir <- args[4]

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

setwd(TT_dir)
temp = list.files(pattern="*.bam$") #Files should be named after the following pattern: TT_TX1072_[Cell-line]_[Day]_[Replicate]_dedup.bam

feature_counts_TT <- featureCounts(temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,
                                GTF.featureType = "gene", strandSpecific = 2, allowMultiOverlap = TRUE)
counts_TT <- data.frame(feature_counts_TT$counts, feature_counts_TT$annotation) %>%
  select(-Chr, -Start, -End, -Strand)
names(counts_TT) <- gsub(x = names(counts_TT), pattern = "\\.", replacement = "_")
names(counts_TT) <- gsub(x = names(counts_TT), pattern = "_dedup_bam", replacement = "")
names(counts_TT) <- gsub(x = names(counts_TT), pattern = "_TX1072", rlseplacement = "")


TPM_TT <- counts_TT %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length) %>% 
  inner_join(gene_names) %>% 
  select(-GeneID) %>%
  select(gene, everything())

setwd(RNA_dir)
temp = list.files(pattern="*.bam$") #Files should be named after the following pattern: RNA_TX1072_[Cell-line]_[Day]_[Replicate]_dedup.bam

feature_counts_RNA <- featureCounts(temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,
                                   GTF.featureType = "exon", strandSpecific = 2, allowMultiOverlap = FALSE, nthreads = 4)
counts_RNA <- data.frame(feature_counts_RNA$counts, feature_counts_RNA$annotation) %>%
  select(-Chr, -Start, -End, -Strand)
names(counts_RNA) <- gsub(x = names(counts_RNA), pattern = "\\.", replacement = "_")
names(counts_RNA) <- gsub(x = names(counts_RNA), pattern = "_dedup_bam", replacement = "")
names(counts_RNA) <- gsub(x = names(counts_RNA), pattern = "_TX1072", replacement = "")

TPM_RNA <- counts_RNA %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length) %>% 
  inner_join(gene_names) %>% 
  select(-GeneID) %>%
  select(gene, everything())

setwd(output_dir)

write_delim(counts_RNA, "rnaseq_counts.txt", delim = "\t")
write_delim(counts_TT, "ttseq_counts.txt", delim = "\t")
write_delim(TPM_RNA, "TPM_rnaseq.txt", delim = "\t")
write_delim(TPM_TT, "TPM_ttseq.txt", delim = "\t")




