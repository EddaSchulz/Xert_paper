#Script to plot Sup Fig. 1b-d
library(tidyverse)
library(UpSetR)
library(gridExtra)
library(egg)
library(extrafont)
loadfonts()

theme_set( theme_classic() +
             theme(legend.title = element_blank(), axis.title = element_text(size = 6, family = "Arial"), legend.position = "none",
                   panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(),
                   plot.title = element_blank(), axis.text = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
char_res <- args[1]
output_dir <- args[2]

char_res_df <- read.delim(char_res)

setwd(output_dir)

upset_df <- char_res_df %>%
  select(id, atac, starr, fantom5) %>%
  mutate_at(vars(-id), funs(ifelse(. == TRUE, 1, 0)))

pdf(file = "Sup_Fig_1b_upsetR_library.pdf", onefile = FALSE, height = unit(3, 'cm'), width = unit(3, 'cm')
    , useDingbats = FALSE)
print(upset(upset_df, order.by = "freq"))
dev.off()

length_plot <- char_res_df  %>%
  ggplot() +
  stat_ecdf(aes(x = length), color = "black", geom = "line", size = 0.5) +
  labs(x = "RE length", y = "Cumulative frequency")


pdf("Sup_Fig_1c_RE_length.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(length_plot, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off() 

nguides_plot <- char_res_df %>%
  ggplot() +
  stat_ecdf(aes(x = nguides), color = "black", geom = "line", size = 0.5) +
  labs(x = "Nr. of sgRNAs / RE", y = "Cumulative frequency") +
  xlim(0, 300)

pdf("Sup_Fig_1d_RE_nguides.pdf", onefile = FALSE, useDingbats = FALSE)
fix <- set_panel_size(nguides_plot, width = unit(2, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()