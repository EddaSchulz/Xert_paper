#Plots embryo rna-seq data from Zhang et al. 2018 for Fig. 4i and Supplementary Fig. 4f
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
            theme(legend.title = element_blank(), axis.title.x = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6, family = "Arial"), axis.title.y = element_text(size = 6, family = "Arial"),
                  plot.title = element_blank(), legend.position = "none", strip.background = element_blank(),
                  strip.text = element_text(size = 6, family = "Arial"), axis.text.x = element_text(angle = 90, vjust = 0.5)))

args <- commandArgs(trailingOnly = TRUE)
TPM_Zhang <- args[1] #Supply path to TPM tables from Zhang et al. 2018 (generated with NGS_downstream/scripts/TPM_Zhang_2018.R)
output_dir <- args[2]

setwd(output_dir)

TPM_embryo <- read.delim(TPM_Zhang)

TPM_long_embryo <- TPM_embryo %>%
  pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
  separate(sample, c("tissue", "replicate"), sep = ".r")


genes <- c("Xert","Pou5f1" ,"Otx2" ,"Ftx" ,"Jpx" ,"Rlim")

tpm_plot <- TPM_long_embryo %>%
  filter(gene %in% genes) %>%
  ggplot(aes(x = tissue, y = tpm, color = factor(gene, levels = genes))) +
  facet_wrap(~factor(gene, levels = genes), scales = "free_y", ncol = 1) +
  geom_point(size = 0.6, alpha = 0.4) +
  stat_summary(fun = "mean", geom = "crossbar", lwd = 0.25, width = 0.5) +
  scale_color_manual(values = c("#00A69C","#AD4994", "#AD4994",  "#EA5B2B", "#D7CB48", "#929497")) +
  labs(y = "TPM") +
  ylim(0, NA)

pdf("Fig_4i_Sup_Fig_4f_TPM_embryo.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(tpm_plot, width = unit(3, "cm"), height = unit(1.25, "cm"))

print(grid.arrange(fix))

dev.off()

