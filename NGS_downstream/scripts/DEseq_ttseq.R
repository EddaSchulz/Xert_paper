#This script calculates DE genes for Supplementary Table 4
library(DESeq2)
library(tidyverse)
library(Rgb)

args <- commandArgs(trailingOnly = TRUE)
TT_counts <- args[1] #Supply path to TT-seq raw counts
RNA_counts <- args[2] #Supply path to RNA-seq raw counts
gencode_xert <- args[3] #Supply path to GTF file containing GENCODE mm10 annotation with Xert added 
output_dir <- args[4]

counts_TT <- read.delim(TT_counts) %>% 
  column_to_rownames("GeneID")
counts_RNA <- read.delim(RNA_counts) %>% 
  column_to_rownames("GeneID")

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()


setwd(output_dir)

# Calculates differential expression between XXdXic and XO at the different timepoints
deseq_fun <- function(counts, day) {
  line <- c("XO", "XO", "XXdXic", "XXdXic")
  timepoint <- rep(day, 4)
  coldata <- data.frame(line, day)
  rownames(coldata) <- colnames(counts)
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ line)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  res <- data.frame(results(DESeq(dds))) %>%
    rownames_to_column("GeneID") %>%
    inner_join(gene_names) %>%
    mutate(day = day)
  return(res)
}

d0_TT_counts <- counts_TT[,c("TT_XO_d0_r1", "TT_XO_d0_r2", "TT_XXdXic_d0_r1", "TT_XXdXic_d0_r2")]
d2_TT_counts <- counts_TT[,c("TT_XO_d2_r1", "TT_XO_d2_r2", "TT_XXdXic_d2_r1", "TT_XXdXic_d2_r2")]
d4_TT_counts <- counts_TT[,c("TT_XO_d4_r1", "TT_XO_d4_r2", "TT_XXdXic_d4_r1", "TT_XXdXic_d4_r2")]

res_d0_TT <- deseq_fun(d0_TT_counts, "d0")
res_d2_TT <- deseq_fun(d2_TT_counts, "d2")
res_d4_TT <- deseq_fun(d4_TT_counts, "d4")
res_total_TT <- rbind(res_d0_TT, res_d2_TT, res_d4_TT) %>%
  mutate(seq = "TT")

d0_RNA_counts <- counts_RNA[,c("RNA_XO_d0_r1", "RNA_XO_d0_r2", "RNA_XXdXic_d0_r1", "RNA_XXdXic_d0_r2")]
d2_RNA_counts <- counts_RNA[,c("RNA_XO_d2_r1", "RNA_XO_d2_r2", "RNA_XXdXic_d2_r1", "RNA_XXdXic_d2_r2")]
d4_RNA_counts <- counts_RNA[,c("RNA_XO_d4_r1", "RNA_XO_d4_r2", "RNA_XXdXic_d4_r1", "RNA_XXdXic_d4_r2")]

res_d0_RNA <- deseq_fun(d0_RNA_counts, "d0")
res_d2_RNA <- deseq_fun(d2_RNA_counts, "d2")
res_d4_RNA <- deseq_fun(d4_RNA_counts, "d4")
res_total_RNA <- rbind(res_d0_RNA, res_d2_RNA, res_d4_RNA) %>%
  mutate(seq = "RNA")

res_total <- rbind(res_total_RNA, res_total_TT)

write_delim(res_total, "DEseq2_ttseq_total.txt", delim = "\t")
