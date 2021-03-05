#This script plots Supplementary Fig. 2g
library(tidyverse)
library(egg)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
  theme(legend.title = element_blank(), 
                panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                axis.text = element_text(size = 6, family = "Arial"), axis.title = element_text(size = 6, family = "Arial"), 
                legend.text = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
count_matrix <- args[1] #Supply path to the count matrix file
output_dir <- args[2]

setwd(output_dir)

#performs PCA analysis and plots first two principal components
count_df <- read.delim(count_matrix) %>%
  select(-c(1:3))

colnames(count_df) <- gsub(".*Zylicz", "nChIP", colnames(count_df))
colnames(count_df) <- gsub(".*CUT", "CUT", colnames(count_df))
colnames(count_df) <- gsub("_TX1072.*H", "_H", colnames(count_df))
colnames(count_df) <- gsub("_d.*", "", colnames(count_df))

pca <- prcomp(count_df, center = TRUE, scale = TRUE)

pca_plot <- data.frame(pca$rotation) %>%
  rownames_to_column(var = "names") %>%
  separate(names, c("method", "mark"), sep = "_") %>%
  ggplot() +
  geom_point(aes(x  = PC1, y = PC2, shape = factor(method, levels = c("CUTnTag", "nChIP")), 
                 color = factor(mark, levels = c("H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3", "H2AK119ub")))) +
  scale_color_manual(values = c("#43AD4F", "#415CA7", "#62C3D4", "#B71398", "#FFD700"))

pdf("Sup_Fig_2g_nChIP_CUTnTag_PCA.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(pca_plot, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()
