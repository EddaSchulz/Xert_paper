#This script plots the correlation matrix for Supplementary Fig. 2e
library(tidyverse)
library(egg)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
  theme(legend.title = element_text(size = 6, family = "Arial"), axis.title = element_blank(),
                panel.border = element_rect(size = 0.5, color = "black", fill = NA), legend.key.size = unit(0.5,"line"),
                axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
                legend.text = element_text(size = 6, family = "Arial"), axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5)))

args <- commandArgs(trailingOnly = TRUE)
input_matrix <- args[1] #Supply path to txt file containing pearson correlation between the merged CUT&Tag samples
out_dir <- args[2]

setwd(out_dir)

pearson_df <- read.delim(input_matrix, skip = 1) %>% 
  rename(sample = X) %>% 
  mutate(sample = gsub(".*_X", "X", sample))  %>% 
  mutate(sample = gsub("_m.*", "", sample))
  
colnames(pearson_df) <- gsub(".*_X", "X", colnames(pearson_df))
colnames(pearson_df) <- gsub("_m.*", "", colnames(pearson_df))

pearson_df_long <- pearson_df %>%
  pivot_longer(-sample, names_to = "sample2", values_to = "pearson")

hclust <- pearson_df %>%
  dplyr::select(-sample)

data <- scale(t(hclust))
ord <- hclust(dist(data, method = "euclidean"), method = "ward.D")$order

sample_levels <- unique(pearson_df$sample)[ord]

pearson_plot <- pearson_df_long %>%
  ggplot(aes(x = factor(sample, levels = sample_levels), y = factor(sample2, levels = sample_levels), fill = pearson)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1), midpoint = 0) +
  labs(fill = "Pearson correlation")

pdf("Sup_Fig_2e_CUTnTag_correlation.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(pearson_plot, width = unit(3, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()
