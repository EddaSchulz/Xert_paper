# Plots Cistrome DB toolkit analysis for Figure S2L
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
            theme(legend.title = element_blank(), axis.title.y = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6, family = "Arial"), axis.title.x = element_text(size = 6, family = "Arial"),
                  plot.title = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(size = 6, family = "Arial"), legend.text = element_text(size = 6, family = "Arial")))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] 
output_dir <- args[2]

setwd(input_dir)

cistrome <- read.csv("Cistrome_toolkit_Xert.csv") 

setwd(output_dir)


wang_2017 <- c("GSM1782914", "GSM1782918")
buecker_2014 <- c("GSM1355169", "GSM1355170", "GSM1355167", "GSM1355168")

cistrome_max <- cistrome %>%
  group_by(Factor) %>%
  summarize(max_GIGGLE = max(GIGGLE_score)) %>%
  slice_max(order_by = max_GIGGLE, n = 15)

cistrome_df <- inner_join(cistrome, cistrome_max) %>%
  mutate(col = ifelse(GSM_ID %in% wang_2017, "Wang 2017", 
                      ifelse(GSM_ID %in% buecker_2014, "Buecker 2014", "-")))

cistrome_plot <- cistrome_df %>%
  ggplot(aes(y = fct_reorder(Factor, GIGGLE_score, max), x = GIGGLE_score, color = col)) +
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("black", "deeppink1", "chartreuse1"), breaks = c("Buecker 2014", "Wang 2017")) +
  labs(x = "GIGGLE score (Similarity)")


pdf("S2L_CistromeDB_Xert.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(cistrome_plot, width = unit(1.5, "cm"), height = unit(3, "cm"))

print(grid.arrange(fix))

dev.off()
