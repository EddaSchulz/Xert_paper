#This script plots the heatmaps for Figure 2J
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
buecker_dir <- args[1]
wang_dir <- args[2]
candidate_re <- args[3]
output_dir <- args[4]
sig_RE = c("RE45", "RE46", "RE47", "RE49", "RE50", "RE51", "RE52", "RE53", "RE55", "RE56", "RE57", "RE58", "RE59", 
           "RE60", "RE61", "RE71", "RE85", "RE87", "RE93", "RE95", "RE96", "RE97")  
           
setwd(output_dir)           

buecker_reads <- list.files(buecker_dir, pattern = ".bam$", full.names = TRUE) %>% 
  str_subset("Input", negate = TRUE)

wang_reads <- list.files(wang_dir, pattern = ".bam$", full.names = TRUE)


feature_counts_buecker <- featureCounts(buecker_reads, annot.ext = candidate_re, isPairedEnd = FALSE)
feature_counts_wang <- featureCounts(wang_reads, annot.ext = candidate_re, isPairedEnd = TRUE)



counts_df <- data.frame(feature_counts_buecker$counts, feature_counts_wang$counts)
names(counts_df) <- gsub(x = names(counts_df), pattern = "\\.", replacement = "_")
names(counts_df) <- gsub(x = names(counts_df), pattern = "2_3", replacement = "2/3")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_dedup_bam", replacement = "")

total_reads <- left_join(feature_counts_buecker$stat, feature_counts_wang$stat) %>%
  dplyr::select(-Status) %>%
  summarize_all(sum) %>%
  mutate_all(funs(. / 1000000))

CPM_df <- sweep(counts_df, 2, as.vector(t(unlist(total_reads))), FUN = '/') %>%
  rownames_to_column("RE") %>%
  filter(RE %in% sig_RE) %>%
  pivot_longer(-RE, names_to = "sample", values_to = "cpm") %>%
  separate(sample, c("tf", "source", "type", "rep"), sep = "_") %>% 
  mutate(type = ifelse(type == "ESC", "undiff", "diff")) %>% 
  group_by(RE, tf, source, type) %>% 
  summarize(cpm = mean(cpm))

lfc <- CPM_df %>%
  pivot_wider(names_from = type, values_from = cpm) %>%
  mutate(lfc = log2(diff/undiff))

write_delim(lfc, "2J_tf_CPM_lfc.txt", delim = "\t")

lfc_map <- lfc %>%
  ggplot(aes(y = factor(tf, levels = c("TCF3", "SMAD2/3", "OCT4", "OTX2")), x = RE)) +
  geom_tile(aes(fill = lfc), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-6, 6)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(fill = "Fold change Diff/Undiff (Log2)")

pdf("2J_tf_lfc_heatmap.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(lfc_map, width = unit(4.587, "cm"), height = unit(0.834, "cm"))

print(grid.arrange(fix))

dev.off()

