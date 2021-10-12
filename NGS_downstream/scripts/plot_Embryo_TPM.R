#Plots embryo rna-seq data in Figures 3I and S3D
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
TPM_table <- args[1] #Supply path to TPM tables from embryo RNA-seq data
output_dir <- args[2]

setwd(output_dir)

TPM_embryo <- read.delim(TPM_table)

TPM_long_embryo <- TPM_embryo %>%
  pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
  mutate(sample = gsub(x = sample, pattern = "5_", replacement = "5-")) %>% 
  mutate(sample = gsub(x = sample, pattern = "0_", replacement = "0-")) %>% 
  separate(sample, c("seq", "source", "tissue", "replicate"), sep = "_")


genes <- c("Xert","Pou5f1" ,"Otx2" ,"Ftx" ,"Jpx" ,"Rlim")
order <- c("Zygote", "2cell", "4cell", "8cell", "16cell", "E3.5-ICM", "E3.5-TE", "E4.0-ICM",
           "E5.5-EPI", "E5.5-VE", "E6.5-EPI", "E6.5-VE", "E7.5-ECT", "E7.5-END", "E7.5-MES",
           "E7.5-PS")

tpm_plot <- TPM_long_embryo %>%
  filter(gene %in% genes) %>%
  ggplot(aes(x = tissue, y = tpm, color = factor(gene, levels = genes))) +
  facet_wrap(~factor(gene, levels = genes), scales = "free_y", ncol = 1) +
  geom_point(size = 0.6, alpha = 0.4) +
  stat_summary(fun = "mean", geom = "crossbar", lwd = 0.25, width = 0.5) +
  scale_color_manual(values = c("#00A69C","#AD4994", "#AD4994",  "#EA5B2B", "#D7CB48", "#929497")) +
  labs(y = "TPM") +
  scale_x_discrete(limits = order) +
  ylim(0, NA)

pdf("3I_S3D_TPM_embryo.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(tpm_plot, width = unit(3.5, "cm"), height = unit(1.25, "cm"))

print(grid.arrange(fix))

dev.off()

