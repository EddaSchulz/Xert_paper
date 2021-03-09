#Plots Sup. Fig. 1h
library(tidyverse)
library(plyr)
library(egg)
library(gridExtra)
library(extrafont)
loadfonts()
summarize <- dplyr::summarize

theme_set(theme_classic()+
            theme(axis.title = element_text(size = 6, family = "Arial"), 
                  panel.border = element_rect(color = "black", size = 0.5, fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"), 
                  strip.background = element_blank(), strip.text = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] 
output_dir <- args[2]

setwd(input_dir)

counts <- read.delim("norm_counts.txt")

setwd(output_dir)

counts_long <- counts %>%
  pivot_longer(c(-id, -re), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("replicate", "fraction")) %>%
  pivot_wider(names_from = replicate, values_from = counts)

counts_mat <- counts %>%
  select(-id, -re) %>%
  as.matrix()

cor_mat <- round(cor(counts_mat),2)

cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sample2", values_to = "cor") %>%
  separate(sample, c("replicate", "fraction"), sep = "_") %>%
  separate(sample2, c("replicate2", "fraction2"), sep = "_") %>%
  filter(replicate != replicate2 & fraction == fraction2) %>%
  select(-replicate, -replicate2, -fraction2) %>%
  unique()

limits <- counts_long %>%
  group_by(fraction) %>%
  summarize(lim = round_any(max(max(R1), max(R2)), 100, ceiling)) 

help_df <- left_join(cor_long, limits)

cor_plot <- counts_long %>%
  ggplot(aes(x = R1, y = R2)) +
  facet_wrap(~factor(fraction, levels = c("Unsorted", "Negative", "Low", "Medium", "High")), scales = "free", nrow = 1) +
  geom_point(size = 0.2, alpha = 0.05) +
  geom_text(data = help_df, aes(x = 0.2 * lim, y = 0.9 * lim, label = paste("r =", cor, sep = " ")), size = 6 / 2.8, 
            family = "Arial") +
  geom_blank(data = help_df, aes(x = lim, y = lim)) +
  lims(y = c(0, NA), x = c(0, NA)) +
  labs(x = "Norm. counts (R1)", y = "Norm. counts (R2)")

pdf("Sup_Fig_1h_screen_correlation.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(cor_plot, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))
dev.off()