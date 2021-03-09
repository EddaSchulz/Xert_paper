# Plots Supplementary Figs. 1i-j
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
    theme(legend.title = element_blank(), legend.text = element_text(size = 6, family = "Arial"),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
          strip.background = element_blank(), axis.title = element_text(size = 6, family = "Arial"),  
          axis.text = element_text(size = 6, family = "Arial"), strip.text = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

setwd(input_dir)


counts <- read.delim("norm_counts.txt")
lib_counts <- read.delim("plasmid_lib_counts.txt")


setwd(output_dir)

counts_long <-counts %>%
  pivot_longer(-c(id, re), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("replicate", "fraction"))

unsorted_counts <- counts_long %>%
  filter(fraction == "Unsorted") %>%
  pivot_wider(names_from = fraction, values_from = counts)

cum_plot_df <- counts_long %>%
  filter(fraction != "Unsorted") %>%
  left_join(unsorted_counts) %>%
  rename(Sorted = counts) %>%
  pivot_longer(c(Sorted, Unsorted), names_to = "curve", values_to = "counts")

cumulative_plot <- cum_plot_df %>%
  ggplot(aes(x = log2(counts), color = curve)) +
  facet_wrap(replicate ~ factor(fraction, levels = c("Negative", "Low", "Medium", "High")), ncol = 4) +
  stat_ecdf(geom = "line", size = 0.3) +
  scale_color_manual(values = c("red", "black")) +
  xlim(5, 11) +
  labs(x = "Normalized sgRNA counts (log2)", y = "Cumulative sgRNA frequency")

pdf("Sup_Fig_1i_cumulative_frequencies.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(cumulative_plot, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

cumulative_plasmid <- lib_counts %>%
  ggplot(aes(x = log2(Plasmid_library))) +
  stat_ecdf(geom = "line", size = 0.3, color = "black") +
  xlim(5, 11) +
  labs(x = "Normalized sgRNA counts (log2)", y = "Cumulative sgRNA frequency")

pdf("Sup_Fig_1e_cumulative_plasmid_library.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(cumulative_plasmid, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

width_df <- left_join(counts, lib_counts) %>%
  filter(re != "NT") %>%
  pivot_longer(c(-id, -re), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)

sample_levels = rev(c("Plasmid_library", "R1_Unsorted", "R1_Negative", "R1_Low", "R1_Medium", "R1_High", 
                  "R2_Unsorted", "R2_Negative", "R2_Low", "R2_Medium", "R2_High"))

width_plot <- width_df %>%
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point(size = 2) +
  labs(x = "Distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1.5, 3)

pdf("Sup_Fig_1j_log2_distribution_width.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(width_plot, width = unit(2, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()

ntc_width_df <- left_join(counts, lib_counts) %>%
  filter(re == "NT") %>%
  pivot_longer(c(-id, -re), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)

sample_levels = rev(c("Plasmid_library", "R1_Unsorted", "R1_Negative", "R1_Low", "R1_Medium", "R1_High", 
                      "R2_Unsorted", "R2_Negative", "R2_Low", "R2_Medium", "R2_High"))

ntc_width_plot <- ntc_width_df %>%
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point(size = 2) +
  labs(x = "NTC distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1, 2.5)


pdf("Sup_Fig_1j_ntc_distribution_width.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(ntc_width_plot, width = unit(2, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()
