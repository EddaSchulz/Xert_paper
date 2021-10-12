#This script performs the bootstrapping analysis
#Plots Figures 1G and S1M
library(tidyverse)
library(gridExtra)
library(egg)
library(matrixStats)
library(extrafont)
loadfonts()

theme_set(
  theme_classic() +
    theme(axis.title = element_blank(), panel.border = element_rect(color = "black", fill=NA, size=0.5),
          axis.line = element_blank(), legend.position = "none",
          axis.ticks = element_blank(), axis.text.y = element_text(margin = margin(r = 0), family = "Arial", size = 6),
          strip.text = element_text(size = 6, family = "Arial"), strip.background = element_blank(),
          axis.text.x = element_blank(), legend.key.size = unit(0.5,"line"))
)

args <- commandArgs(trailingOnly = TRUE)
raw_counts <- args[1] #Counts are found at /CRISPR_screen_analysis/files/raw_counts.txt
output_dir <- args[2]

setwd(output_dir)

counts <- read.delim(raw_counts) %>%
  filter(R1_Unsorted > 10 | R2_Unsorted > 10)

counts_long <- counts %>% pivot_longer(c(-id, -re), names_to = "bin", values_to = "counts")
nt_counts <- counts %>% filter(re == "NT") %>%
  group_by(re) %>% summarise_if(.predicate = function(x) is.numeric(x), .funs = funs(mean)) %>%
  pivot_longer(c(-re), names_to = "bin", values_to = "NT_counts") %>% select(-re)
counts_norm <- merge(counts_long, nt_counts) %>%
  mutate(counts = counts / NT_counts * 1000) %>%
  select(-NT_counts) %>%
  pivot_wider(names_from = "bin", values_from = "counts")

nr_sg <- counts %>%
  group_by(re) %>%
  summarize(nr_sg = length(id))

mean_counts <- counts_norm %>%
  pivot_longer(c(-id, -re), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("rep", "fraction")) %>%
  group_by(id, re, fraction) %>%
  summarize(mean_counts = mean(counts))

rel_counts <- mean_counts %>%
  pivot_wider(names_from = "fraction", values_from = "mean_counts") %>%
  mutate_at(vars(Negative, Low, Medium, High), funs(. / Unsorted)) %>%
  select(id, re, Negative, Low, Medium, High, -Unsorted)

boot_size = 50
nr_cyc = 1000
boot_all=data.frame()

sel_re = nr_sg %>%
  filter(nr_sg > boot_size)


for (c in 1:length(sel_re$re)) {
  t <- sel_re$re[c]
  sel_data <- rel_counts %>%
    filter(re == t)
  nr_guides <- dim(sel_data)[1]

  for (i in 1:nr_cyc) {
    sel_sg <- ceiling(runif(boot_size, min=0, max=nr_guides))
    sample_mean <- colMeans(sel_data[sel_sg,3:6])
    if (i==1) {
      boot_data <- sample_mean
    } else {
      boot_data <- rbind(boot_data, sample_mean)
    }
  }
  boot_long = data.frame(it = 1:nr_cyc, boot_data) %>%
    pivot_longer(-it, names_to = "fraction", values_to = "fc") %>%
    mutate(re = t)
  boot_all = rbind(boot_all, boot_long)
}



empirical_pvalue <- boot_all %>%
  mutate(n_high = ifelse(fc > 1, 1, 0), n_low = ifelse(fc < 1, 1, 0)) %>%
  group_by(fraction, re) %>%
  summarize(p_high = (1000 - sum(n_high)) / 1000, p_low = (1000 - sum(n_low)) / 1000,
            num_high = sum(n_high), num_low = sum(n_low)) %>%
  mutate(padj_high = p.adjust(p_high, "fdr", n = nrow(sel_re)), padj_low = p.adjust(p_low, "fdr", n = nrow(sel_re)))

sig_fdr <- empirical_pvalue %>%
  mutate(p_min = pmin(padj_high, padj_low)) %>%
  mutate(sig = ifelse(p_min <= 0.01, TRUE, FALSE))

boot_out = left_join(boot_all, sig_fdr) %>%
  mutate(log_fc = log2(fc))



help_df <- boot_out %>%
  group_by(re) %>%
  summarize(max = max(log_fc), min = min(log_fc)) %>%
  mutate(lim = 0.1 * ceiling(10 * pmax(max * 1.25, -(0.9 * min)))) %>%
  mutate(y_pos = 0.88 * lim)

plot_df <- inner_join(boot_out, help_df) %>%
  mutate(fraction = factor(fraction, levels = c("Negative", "Low", "Medium", "High")))
text_df <- plot_df %>%
  group_by(fraction, re, sig, y_pos) %>%
  summarize()

boot_plot <- plot_df %>%
  ggplot() +
  facet_wrap(~re, scales = "free", ncol = 10) +
  geom_hline(yintercept = 0, linetype = "dotted",size = 0.5) +
  geom_boxplot(aes(x = fraction, y = log_fc, color = fraction), outlier.size = 0.1, lwd = 0.1) +
  geom_boxplot(aes(x = fraction, y = log_fc, fill = fraction), outlier.shape = NA, lwd = 0.1) +
  geom_text(data = subset(text_df, sig == TRUE), aes(x = fraction , y = y_pos , label='*'), size=4, color='black')  +
  geom_blank(aes(x = fraction, y = lim)) +
  geom_blank(aes(x = fraction, y = -lim)) +
  scale_fill_grey(start = 0.3, end = 0.9) +
  scale_color_grey(start = 0.3, end = 0.9)

pdf("1G_S1M_bootstrapping_plot.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(boot_plot, width = unit(1, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()

pval_df <- boot_out %>%
  select(fraction, re, p_high, p_low, padj_high, padj_low) %>%
  unique()

write_delim(pval_df, "bootstrapping_pval.txt", delim = "\t")

rank_df <- boot_out %>%
  group_by(fraction, re) %>%
  summarize(log_fc = mean(log_fc)) %>%
  pivot_wider(names_from = fraction, values_from = log_fc) %>%
  mutate(mean_fc = (High + Medium + Low - Negative) /  4) %>%
  arrange(mean_fc)

write_delim(rank_df, "bootstrapping_rank.txt", delim = "\t")
