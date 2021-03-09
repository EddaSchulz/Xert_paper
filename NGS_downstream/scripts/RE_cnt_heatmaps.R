#This script plots the heatmaps in figures 3C and 3D
library(tidyverse)
library(Rsubread)
library(extrafont)
library(gridExtra)
library(EnvStats)
library(egg)
loadfonts()

theme_set(theme_classic() +
  theme(legend.title = element_text(size = 6, family = "Arial"), plot.title = element_text(hjust = 0.5),
                axis.title = element_blank(), panel.border = element_rect(size = 0.5, color = "black", fill = NA),
                axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
                legend.text = element_text(size = 6, family = "Arial"), axis.ticks = element_blank(),
                strip.background = element_blank(), strip.text = element_text(size = 6, family = "Arial"),
                legend.key.size = unit(0.5,"line")))

args <- commandArgs(trailingOnly = TRUE)
atac_dir <- args[1] #Supply path to directory containing ATAC-seq BAM files (generated via /NGS_alignment/ATACseq_align.sh)
cnt_dir <- args[2] #Supply path to directory containing CUT&tag BAM files (generated via /NGS_alignment/CUTnTag_align.sh)
candidate_re <- args[3] #Supply path to SAF file containing candidate RE's (At /NGS_downstream/files/candidate_re.saf)
output_dir <- args[4] #Supply path to output file
sig_RE = c("RE58", "RE57", "RE96", "RE59", "RE85", "RE93", "RE61", "RE97", "RE95", "RE60", "RE71", "RE87",
           "RE56", "RE22", "RE55", "RE53", "RE12", "RE45", "RE47", "RE46", "RE50", "RE49", "RE52", "RE51") #List of RE's within chrX:103182701-103955531 that were significant in the MAGeCK analysis (FDR <= 0.05)

ATAC_reads <- list.files(atac_dir, pattern = ".bam$", full.names = TRUE)

CnT_reads <- list.files(cnt_dir, pattern = ".bam$", full.names = TRUE)
CnT_filter <- CnT_reads[str_detect(CnT_reads, "H3K4me1|H3K4me3|H3K27ac")]

cd(output_dir)

temp = c(ATAC_reads, CnT_filter)

feature_counts <- featureCounts(temp, annot.ext = candidate_re, isPairedEnd = TRUE, nthreads = 4)

counts_df <- data.frame(feature_counts$counts)
names(counts_df) <- gsub(x = names(counts_df), pattern = "\\.", replacement = "_")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_dedup_bam", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_sorted_blacklisted_bam", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_TX1072", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "CUTnTag_", replacement = "")

total_reads <- feature_counts$stat %>%
  dplyr::select(-Status) %>%
  summarize_all(sum) %>%
  mutate_all(funs(. / 1000000))

CPM_df <- sweep(counts_df, 2, as.vector(t(unlist(total_reads))), FUN = '/') %>%
  rownames_to_column("RE") %>%
  filter(RE %in% sig_RE) %>%
  pivot_longer(-RE, names_to = "sample", values_to = "cpm") %>%
  separate(sample, c("line", "mark", "day", "replicate"), sep = "_") %>% 
  mutate(line_extra = ifelse(line == "ATAC",mark ,line), mark_extra = ifelse(line =="ATAC", line, mark),
         line = line_extra, mark = mark_extra) %>%
  select(-line_extra, -mark_extra)

log_df <- CPM_df %>%
  group_by(RE, line, mark, day) %>%
  summarize(cpm = geoMean(cpm + 0.01)) %>%
  mutate(log_cpm = log2(cpm))

help_df <- log_df %>%
  group_by(RE, mark) %>%
  summarize(mean_cpm = mean(log_cpm), sd_cpm = sd(log_cpm))

zscore_df <- left_join(log_df, help_df) %>%
  mutate(zscore = (log_cpm - mean_cpm)/sd_cpm)

abs_counts_z <- counts_df %>%
  rownames_to_column("RE") %>%
  filter(RE %in% sig_RE) %>%
  pivot_longer(-RE, names_to = "sample", values_to = "counts") %>%
  separate(sample, c("line", "mark", "day", "replicate"), sep = "_") %>%
  group_by(RE, line, mark, day) %>%
  summarize(counts = min(counts)) %>%
  mutate(sig = ifelse(counts >= 5, 1, 0)) %>%
  group_by(RE, mark) %>%
  summarize(sig = sum(sig) == 0)

zscore_map <- zscore_df %>%
  ggplot(aes(x = factor(mark, levels = c("ATAC", "H3K4me3", "H3K27ac", "H3K4me1")), y = reorder(RE, desc(RE)))) +
  facet_wrap(day~factor(line, levels = c("XXdXic", "XO")), nrow = 1) +
  geom_tile(aes(fill = zscore), color = "black") +
  geom_tile(data = subset(abs_counts_z, sig), fill = "darkgrey", color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-2, 2)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(fill = "Z-score")

pdf("Fig_3c_RE_zscore.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(zscore_map, width = unit(0.8333, "cm"), height = unit(5, "cm"))

print(grid.arrange(fix))

dev.off()

abs_counts_lfc <- counts_df %>%
  rownames_to_column("RE") %>%
  filter(RE %in% sig_RE) %>%
  pivot_longer(-RE, names_to = "sample", values_to = "counts") %>%
  separate(sample, c("line", "mark", "day", "replicate"), sep = "_") %>% 
  mutate(line_extra = ifelse(line == "ATAC",mark ,line), mark_extra = ifelse(line =="ATAC", line, mark),
         line = line_extra, mark = mark_extra) %>%
  select(-line_extra, -mark_extra) %>% 
  group_by(RE, line, mark, day) %>%
  summarize(counts = min(counts)) %>%
  mutate(sig = ifelse(counts >= 5, 1, 0)) %>%
  group_by(RE, mark, day) %>%
  summarize(sig = sum(sig) == 0)

lfc <- left_join(CPM_df, abs_counts_lfc) %>%
  group_by(RE, line, mark, day, sig) %>%
  summarize(cpm = geoMean(cpm + 0.01)) %>%
  pivot_wider(names_from = line, values_from = cpm) %>%
  mutate(lfc = log2(XXdXic/XO))

lfc_map <- lfc %>%
  ggplot(aes(x = factor(mark, levels = c("ATAC", "H3K4me3", "H3K27ac", "H3K4me1")), y = reorder(RE, desc(RE)))) +
  facet_wrap(~day, nrow = 1) +
  geom_tile(aes(fill = lfc), color = "black") +
  geom_tile(data = subset(lfc, sig), fill = "darkgrey", color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-7, 7)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(fill = "Fold change XXdXic / XO (Log2)")

pdf("Fig_3d_RE_LFC.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(lfc_map, width = unit(0.8333, "cm"), height = unit(5, "cm"))

print(grid.arrange(fix))

dev.off()

