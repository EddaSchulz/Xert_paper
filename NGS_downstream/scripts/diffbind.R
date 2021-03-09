#Runs the Diffbind analysis for differential peak calling of the CUT&Tag and ATAC-seq data
library(tidyverse)
library(DiffBind)

args <- commandArgs(trailingOnly = TRUE)
peak_dir <- args[1] #Supply path to directory containing noX peak files
atac_dir <- args[2] #Supply path to directory containing BAM files (ATAC-seq)
cnt_dir <- args[3] #Supply path to directory containing BAM files (CUT&Tag)
output_dir <- args[4]

#Read in peak file paths and store in a dataframe
Peaks <- list.files(peak_dir, full.names = TRUE)
peak_df <- data.frame(Peaks) %>%
  mutate(SampleID = str_extract(Peaks, "[A-Za-z0-9_]*(?=_noX\\.bed)"))

#Read in BAM file paths and combine with peak file dataframe
bamReads <- list.files(atac_dir, pattern = ".bam$", full.names = TRUE)
atac_prep <- data.frame(bamReads) %>%
  mutate(SampleID = str_extract(bamReads, "[A-Za-z0-9_]*(?=\\_dedup.bam)")) 

bamReads <- list.files(cnt_dir, pattern = ".bam$", full.names = TRUE)
cnt_prep <- data.frame(bamReads) %>%
  mutate(SampleID = str_extract(bamReads, "[A-Za-z0-9_]*(?=\\_sorted_blacklisted.bam)")) %>% 
  filter(str_detect(SampleID, "H3K4me1|H3K4me3|H3K27ac"))

diffbind_prep <- rbind(atac_prep, cnt_prep) %>% 
  inner_join(peak_df) %>%
  mutate(SampleID = gsub("CUTnTag_", "", gsub("TX1072_", "", SampleID))) %>% 
  separate(SampleID, c("Tissue", "Factor", "Condition", "Replicate"), sep = "_") %>% 
  mutate(line_extra = ifelse(Tissue == "ATAC",Factor ,Tissue), mark_extra = ifelse(Tissue =="ATAC", Tissue, Factor),
         Tissue = line_extra, Factor = mark_extra) %>%
  mutate(SampleID = paste(Tissue, Factor, Condition, Replicate, sep = "_"), 
         Treatment = paste(Tissue, Condition, sep = "_")) %>% 
  select(-line_extra, -mark_extra) %>%
  mutate(Treatment = paste(Tissue, Condition, sep = "_")) %>%
  select(SampleID, Tissue, Factor, Condition, Replicate, Treatment, bamReads, Peaks)

#Create dataframe containing all comparisons
comps <- c("clone_d0", "clone_d2", "clone_d4", "diff_d2", "diff_d4")
treatA <- c("XXdXic_d0", "XXdXic_d2", "XXdXic_d4", "XXdXic_d0", "XXdXic_d0")
treatB <- c("XO_d0", "XO_d2", "XO_d4", "XXdXic_d4", "XXdXic_d4")
Factor <- c("ATAC", "H3K4me3", "H3K27ac", "H3K4me1")

comp_treatments <- data.frame(comps, treatA, treatB) 
comp_df <- list(comps = comps, Factor = Factor) %>%
  expand.grid() %>%
  mutate(id = paste(comps, Factor, sep = "_")) %>%
  left_join(comp_treatments)

#Create a list of dataframes containing the necessary data for all comparisons
f_subset <- function(x) {
  subset(diffbind_prep, Factor == x[2] & Treatment %in% c(x[4], x[5]))
}

diffbind_list <- apply(comp_df, 1, f_subset)
names(diffbind_list) <- comp_df$id

#Create DBA objects
diffbind_sheet <- lapply(diffbind_list, function(x) dba(sampleSheet = x))

#Count reads to create a binding matrix
diffbind_count <- lapply(diffbind_sheet, function(x) dba.count(x, bParallel=FALSE))

#Create contrast between conditions
diffbind_contrast <- lapply(diffbind_count, function(x) dba.contrast(x, categories=DBA_TREATMENT, minMembers = 2, bNot = FALSE))

#Perform differential binding analysis
diffbind_analysis <- lapply(diffbind_contrast, function(x) dba.analyze(x, DBA_ALL_METHODS))

#Export dba objects as .txt files
diffbind_tracks <- lapply(diffbind_analysis, function(x) dba.report(x, method=DBA_ALL_METHODS, contrast = 1, th=0.05, bUsePval = F, DBA_DATA_FRAME))
diffbind_df <- lapply(diffbind_tracks, function(x) as.data.frame(x))
sapply(names(diffbind_df), function (x) write_delim(diffbind_df[[x]], paste0(output_dir, x,"_diff.bed"), delim = "\t", col_names = FALSE))

#Export consensus peaks as .txt files
consensus <- lapply(diffbind_sheet, function(x) dba.peakset(x, consensus = DBA_CONDITION, minOverlap = 4, bRetrieve = TRUE))
consensus_df <- lapply(consensus, function(x) as.data.frame(x[,c(1:3)]))
sapply(names(consensus), function (x) write_delim(consensus_df[[x]], paste0(output_dir, x, "_consensus.bed"), delim = "\t", col_names = FALSE))
