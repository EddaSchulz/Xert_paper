#This script plots the emission states of the ChromHMM analysis for Figure S2H
#Requires the emmision state .txt file from /NGS_downstream/ChromHMM.sh as input
library(tidyverse)
library(egg)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
  theme(legend.title = element_text(size = 6, family = "Arial"), plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 6, family = "Arial"), panel.border = element_rect(size = 0.5, color = "black", fill = NA),
          axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
          legend.position = "top", legend.text = element_text(size = 6, family = "Arial"),
          axis.ticks = element_blank(), legend.key.size = unit(0.5,"line"),
          axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank()))

args <- commandArgs(trailingOnly = TRUE)
emission_states <- args[1] #Supply path to .txt file containing emission state signals (example file at NGS_downstream/files/emissions_12.txt)
output_dir <- args[2] #Supply output path
ord <- args[3]

setwd(output_dir)

state_order <- as.integer(unlist(str_split(ord, ",")))

emissions <- read.delim(emission_states) %>%
  rename(state = state..Emission.order.) %>%
  pivot_longer(-state, names_to = "mark", values_to = "signal_strength")

mark_levels = c("ATAC", "H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3")

state_plot <- emissions %>%
  ggplot(aes(x = factor(mark, levels = mark_levels),
             y = factor(state, levels = state_order), fill = signal_strength)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(y = "State (Emission order)", fill = "Signal strength [a.u.]")

pdf("S2H_ChromHMM_States.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(state_plot, width = unit(1.25, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()
