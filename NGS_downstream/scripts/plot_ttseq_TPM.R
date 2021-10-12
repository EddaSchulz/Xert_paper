#This script plots Figure 3B
library(extrafont)
library(gridExtra)
library(egg)
library(tidyverse)
loadfonts()

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6, family = "Arial"),
                  axis.title = element_blank(), panel.border = element_rect(color = "black", size = 0.5, fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
                  legend.text = element_text(size = 6, family = "Arial"), strip.background = element_blank(),
                  strip.text = element_text(size = 6, family = "Arial"), legend.key.size = unit(0.5,"line")))

args <- commandArgs(trailingOnly = TRUE)
TPM_TT_file <- args[1] #Supply path to TT-seq TPM
TPM_RNA_file <- args[2] #Supply path to RNA-seq TPM
DEseq_file <- args[3] #Supply path to total DEseq sheet
output_dir <- args[4]

TPM_TT <- read.delim(TPM_TT_file)
TPM_RNA <- read.delim(TPM_RNA_file)
DEseq_total <- read.delim(DEseq_file)

setwd(output_dir)

TPM_total <- inner_join(TPM_TT, TPM_RNA) %>%
  select(gene, everything())

TPM_long <- TPM_total %>%
  pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
  separate(sample, c("seq", "line", "day", "replicate"), sep = "_")

lnc_genes <- c("Tsix", "Jpx", "Ftx", "Xert", "Xist")

TPM_lnc <- TPM_long %>%
  filter(gene %in% lnc_genes) %>%
  mutate(log_tpm = log2(tpm + 0.1))

#Plots dotplots for the TPM of lncRNA genes within the Xic
TPM_DEseq_lnc <- left_join(TPM_lnc, DEseq_total)
  

#Calculates the position of the significance asterisk
star_df <- TPM_DEseq_lnc %>%
  group_by(gene, seq) %>% 
  summarize(star_pos = max(log_tpm) + 0.075 * (max(log_tpm) - min(log_tpm)), 
            top_border = max(log_tpm) + 0.2 * (max(log_tpm) - min(log_tpm)))

help_df <- TPM_DEseq_lnc %>%
  select(gene, seq, day, padj) %>%
  unique() %>%
  inner_join(star_df) 


TPM_plots <- TPM_DEseq_lnc %>%
  ggplot(aes(x = day)) +
  facet_wrap(seq ~ gene, scales = "free_y", nrow = 2) +
  geom_point(aes(y = log_tpm, color = line), size = 0.5) +
  stat_summary(aes(y = log_tpm, group = line, color = line), 
               geom = "line", fun = "mean", lwd = 0.25) +
  geom_text(data = subset(help_df, padj <= 0.05),
            aes(label = "*", y = star_pos), family = "Arial", color = "black") +
  scale_color_manual(values = c("#5CB4AC", "#EF7E1A")) +
  geom_blank(data = help_df, aes(y = top_border)) +
  labs(y = "TPM + 0.1 (Log2)") +
  scale_y_continuous(breaks = seq(-10, 10, by = 1)) +
  scale_x_discrete(labels = c(0, 2, 4))

pdf("3B_lncRNA_lineplots.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(TPM_plots, height = unit(1.5, "cm"), width = unit(1, "cm"))

print(grid.arrange(fix))

dev.off()


