#Plots Fig. 1d-e and Sup. Fig. 1k-l
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
library(plyr)
loadfonts()

theme_set(theme_classic() +
  theme(legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),  axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 6, family = "Arial"), axis.title.y = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] #Should contain Char_candidate_RE.txt, mageck_mle.txt, mageck_test_RE_high.txt, mageck_test_RE_negative.txt, raw_counts.txt, Char_guides.txt, mageck_test_high.txt, mageck_test_negative.txt   
output_dir <- args[2]
  
setwd(input_dir)

char_RE <- read.delim("Char_candidate_RE.txt") %>% #Can be found at /CRISPR_library_design/files/
  select(re, chr, start, end)
mle <- read.delim("mageck_mle.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  select(re, Negative.beta, Negative.fdr = Negative.wald.fdr, Low.beta, Low.fdr = Low.wald.fdr, 
                Medium.beta, Medium.fdr = Medium.wald.fdr, High.beta, High.fdr = High.wald.fdr)
high_test_re <- read.delim("mageck_test_RE_high.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  select(re, High_LFC = neg.lfc)
negative_test_re <- read.delim("mageck_test_RE_negative.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  select(re, Negative_LFC = neg.lfc)

filter_guides <- read.delim("raw_counts.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  filter(R1_Unsorted >= 10 & R2_Unsorted >= 10 & re != "NT") %>%
  select(id)
char_guides <- read.delim("Char_guides.txt") %>% #Can be found at /CRISPR_library_design/files/
  select(id, re, start, end)
high_test <- read.delim("mageck_test_high.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  select(id, re, p_high = FDR, high_lfc = LFC)
negative_test <- read.delim("mageck_test_negative.txt") %>% #Can be found at /CRISPR_screen_analysis/files/
  select(id, re, p_negative = p.twosided, negative_lfc = LFC)

setwd(output_dir)

re_df <- join_all(list(char_RE, mle, high_test_re, negative_test_re)) %>%
  arrange(desc(High.fdr))
guide_df <- join_all(list(filter_guides, char_guides, high_test, negative_test)) %>%
  arrange(desc(p_high))
fill_high <- ifelse(re_df$High.fdr <= 0.05 & re_df$High_LFC <= 0, "red",
              ifelse(re_df$High.fdr <= 0.05 & re_df$High_LFC >= 0, "blue", "grey"))
color_guides_high <- ifelse(guide_df$p_high <= 0.05 & guide_df$high_lfc <= 0, "red",
                     ifelse(guide_df$p_high <= 0.05 & guide_df$high_lfc >= 0, "blue", "grey"))


xic_plot_high <- ggplot() +
  geom_point(data = guide_df, aes(x = start, y = high_lfc), color = color_guides_high, size = 0.3) +
  geom_abline(aes(slope = 0, intercept = 0)) +
  geom_point(data = re_df, aes(x = (start + end) / 2, y = High_LFC), fill = fill_high, shape = 21, size = 1) +
  labs(y = "Fold change High/Unsorted (log2)") + 
  ylim(-6, 2) +
  scale_x_continuous(limits = c(103198658, 104058961), expand = c(0.005, 0)) 

pdf("Fig_1d_FC_xic_high.pdf",  useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(xic_plot_high, width = unit(16, "cm"), height = unit(5, "cm"))

print(grid.arrange(fix))

dev.off()

xist_plot_high <- re_df %>%
  ggplot() +
  geom_abline(aes(slope = 0, intercept = 1.30103), linetype = "dotted") +
  geom_point(aes(x = (start + end) / 2, y = -log10(High.fdr)), fill = fill_high, shape = 21, size = 1.5) +
  scale_fill_manual(values = c("red", "grey")) +
  labs(y = "FDR High/Unsorted (-log10)")+ 
  scale_x_continuous(limits = c(103430517, 103507425), expand = c(0, 0)) +
  ylim(0, 30)

pdf("Fig_1e_FDR_xist_zoomin_high.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(xist_plot_high, width = unit(6, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()


xert_plot_high <- re_df %>%
  ggplot() +
  geom_abline(aes(slope = 0, intercept = 1.30103), linetype = "dotted") +
  geom_point(aes(x = (start + end) / 2, y = -log10(High.fdr)), fill = fill_high, shape = 21, size = 1.5) +
  scale_fill_manual(values = c("red", "grey")) +
  labs(y = "FDR High/Unsorted (-log10)")+ 
  scale_x_continuous(limits = c(103559910, 103698414), expand = c(0, 0)) +
  ylim(0, 30)

pdf("Fig_1e_FDR_xert_zoomin_high.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(xert_plot_high, width = unit(6, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()


re_sort <- re_df %>%
  arrange(desc(Negative.fdr))
guide_sort <- guide_df %>%
  arrange(desc(p_negative))
fill_negative <- ifelse(re_sort$Negative.fdr <= 0.05 & re_sort$Negative_LFC >= 0, "red",
                    ifelse(re_sort$Negative.fdr <= 0.05 & re_sort$Negative_LFC <= 0, "blue", "grey"))
color_guides_negative <- ifelse(guide_sort$p_negative <= 0.05 & guide_sort$negative_lfc >= 0, "red",
                            ifelse(guide_sort$p_negative <= 0.05 & guide_sort$negative_lfc <= 0, "blue", "grey"))


xic_plot_negative <- ggplot() +
  geom_point(data = guide_sort, aes(x = start, y = negative_lfc), color = color_guides_negative, size = 0.3) +
  geom_abline(aes(slope = 0, intercept = 0)) +
  geom_point(data = re_sort, aes(x = (start + end) / 2, y = Negative_LFC), fill = fill_negative, shape = 21, size = 1) +
  labs(y = "Fold change Negative/Unsorted (log2)") + 
  ylim(-2, 3) +
  scale_x_continuous(limits = c(103198658, 104058961), expand = c(0.005, 0)) 

pdf("Sup_Fig_1k_FC_xic_negative.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(xic_plot_negative, width = unit(16, "cm"), height = unit(5, "cm"))

print(grid.arrange(fix))

dev.off()

xist_plot_negative <- re_sort %>%
  ggplot() +
  geom_abline(aes(slope = 0, intercept = 1.30103), linetype = "dotted") +
  geom_point(aes(x = (start + end) / 2, y = -log10(Negative.fdr)), fill = fill_negative, shape = 21, size = 1.5) +
  scale_fill_manual(values = c("red", "grey")) +
  labs(y = "FDR Negative/Unsorted (-log10)")+ 
  scale_x_continuous(limits = c(103430517, 103507425), expand = c(0, 0)) +
  ylim(0, 35)

pdf("Sup_Fig_1l_FDR_xist_zoomin_negative.pdf", useDingbats  = FALSE, onefile =  FALSE)
fix <- set_panel_size(xist_plot_negative, width = unit(6, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()
