#This script is used to plot the karyotyping heatmaps in Sup. Fig. 5f
library(tidyverse)
library(egg)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6, family = "Arial"), plot.title = element_text(hjust = 0.5),
                  axis.title = element_blank(), panel.border = element_rect(size = 0.5, color = "black", fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
                  legend.text = element_text(size = 6, family = "Arial"), axis.ticks = element_blank(),
                  strip.background = element_blank(), strip.text = element_blank(),
                  legend.key.size = unit(0.5,"line"), axis.text.x = element_text(angle = 90, vjust = 0.5)))
          
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] 
output_dir <- args[2]

setwd(input_dir)

counts <- read.delim("karyotype_counts.txt")


setwd(output_dir)

div_fun <- function(a,b) {
  a / b
}
sum_fun <- function(a) {
  a / sum(a)
}

#The TX1072_XX_A3 cell line is used as the wildtype control
karyo_df <- counts %>%
  select(-start, -end) %>% 
  mutate_at(vars(-chr), ~sum_fun(.)) %>%
  mutate_at(vars(-chr, -TX1072_XX_A3), ~div_fun(., TX1072_XX_A3))  %>%
  mutate(chr = gsub("chr", "", chr)) %>% 
  select(-TX1072_XX_A3) %>%
  pivot_longer(-chr, names_to = "line", values_to = "normed")


chr_order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "X")

karyo_plot <- karyo_df %>%
  ggplot(aes(x = factor(chr, levels = chr_order), y = line, fill = normed)) +
  facet_wrap(~line, scales = "free_y", ncol = 1) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1, limits = c(0.5, 1.5)) +
  labs(fill = "Counts ratio (Clone / XX control)")

pdf("karyo_plot.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(karyo_plot, width = unit(4, "cm"), height = unit(0.2, "cm"))

print(grid.arrange(fix))

dev.off()
