#Script ranks RE's based on beta-score and plots the heatmap in Figure 1F
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6, family = "Arial"),
                  axis.title = element_blank(), panel.border = element_rect(size = 0.5, color = "black", fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"), 
                  axis.text.x = element_text(angle = 90, vjust = 0.5), legend.key.size = unit(0.5,"line"),
                  legend.text = element_text(size = 6, family = "Arial"), axis.ticks = element_blank()))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] #Can be found at /CRISPR_screen_analysis/files/mageck_mle.txt (generated with /CRISPR_screen_analysis/MAGeCK_mle.sh)
output_dir <- args[2]

setwd(input_dir)

mle <- read.delim("mageck_mle.txt") %>%
  select(id = re, Negative.beta, Negative.fdr = Negative.wald.fdr, Low.beta, Low.fdr = Low.wald.fdr, 
         Medium.beta, Medium.fdr = Medium.wald.fdr, High.beta, High.fdr = High.wald.fdr)

setwd(output_dir)

heatmap_df <- mle %>%
  filter(Negative.fdr <= 0.05 | Low.fdr <= 0.05 | Medium.fdr <= 0.05 | High.fdr <= 0.05) %>%
  mutate(score = (-Negative.beta + Low.beta + Medium.beta + High.beta) / 4) %>%
  pivot_longer(cols = c(Negative.beta, Low.beta, Medium.beta, High.beta), names_to = "fraction", values_to = "beta_score") %>%
  pivot_longer(cols = c(Negative.fdr, Low.fdr, Medium.fdr, High.fdr), names_to = "fraction2", values_to = "fdr") %>%
  mutate_at(vars(fraction, fraction2), funs(gsub("\\.[a-z]*", "", .))) %>% 
  filter(fraction == fraction2) %>%
  dplyr::select(-fraction2) %>%
  arrange(score)

heatmap_plot <- heatmap_df %>%
  ggplot(aes(x = factor(fraction, levels = c("Negative", "Low", "Medium", "High")), 
             y = factor(id, levels=unique(id[order(-score)])))) +
  geom_tile(aes(fill = beta_score), color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
  geom_text(data = subset(heatmap_df, fdr <= 0.05), aes(label = "*"), family = "Arial", size = 2.14, nudge_y = -0.15) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Beta score") 

pdf("1F_betascore_heatmap.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(heatmap_plot, width = unit(0.92, "cm"), height = unit(6, "cm"))

print(grid.arrange(fix))

dev.off()
