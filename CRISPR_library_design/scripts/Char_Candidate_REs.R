#Script generates Char_candidate_RE.txt file for further use
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

setwd(input_dir)

candidate_res <- import("Candidate_RE.bed") #Can be found at /CRISPR_library_design/files/
atac <- import("ATAC_peaks_XIC.bed") #Can be found at /CRISPR_library_design/files/ (generated with /CRISPR_library_design/ATAC_peaks.sh)
starr <- import("STARR_peaks_XIC.bed") #Can be found at /CRISPR_library_design/files/ (generated with /CRISPR_library_design/STARR_peaks.sh)
fantom5 <- import("F5_enhancers_XIC.bed") #Can be found at /CRISPR_library_design/files/ (generated with /CRISPR_library_design/FANTOM5_XIC.sh)
guides <- import("gRNAs.bed") #Can be found at /CRISPR_library_design/files/

setwd(output_dir)

atac_vec <- findOverlaps(candidate_res, atac, select = "first")
starr_vec <- findOverlaps(candidate_res, starr, select = "first")
fantom5_vec <- findOverlaps(candidate_res, fantom5, select = "first")

atac_col <- !is.na(atac_vec)
starr_col <- !is.na(starr_vec)
fantom5_col <- !is.na(fantom5_vec)
guide_col <- countOverlaps(candidate_res, guides)
id_col <- paste("RE", c(1:138), sep = "")

char_candidate_res <- read.delim(paste0(input_dir, "/Candidate_RE.bed"), header = FALSE) %>%
  transmute(id = id_col, chr = V1, start = V2, end = V3, atac = atac_col, starr = starr_col, fantom5 = fantom5_col,
            length = V3 - V2, nguides = guide_col)

write_delim(char_candidate_res, "Char_candidate_RE.txt", delim = "\t")
