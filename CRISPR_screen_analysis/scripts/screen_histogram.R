#Script analyzes the FlowFISH data from the screen and plots histograms for figure S1E
library(flowCore)
library(openCyto)
library(tidyverse)
library(egg)
library(gridExtra)
library(extrafont)
library(scales)
loadfonts()
filter <- dplyr::filter

theme_set(theme_classic() + 
    theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial"),
          axis.title = element_text(size = 6, family = "Arial"), strip.text = element_text(size = 6, family = "Arial"),
          strip.background = element_blank()))

args <- commandArgs(trailingOnly = TRUE)
fcs_dir <- args[1] 
output_dir <- args[2] 

flowset <- read.flowSet(path = fcs_dir)
setwd(output_dir)

rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.1e5, 2e5), "SSC-A"= c(0.3e5, 2e5))
BoundaryFrame <- Subset (flowset[[1]], rectGate) 

sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"))
singlets_flowset <- Subset (flowset, rectGate %&% sing_g)

extr_frame <- fsApply(flowset, exprs, simplify = FALSE)
extr_df <- lapply(extr_frame, data.frame)
flow_df <- data.frame(do.call(rbind,extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub(".fcs.*", "", sample), Xist = APC...670.14...A) %>%
  filter(Xist > 0) %>%
  separate(sample, c("replicate", "sample"), sep = "_")

facs_histogram <- ggplot() +
  facet_wrap(~replicate) +
  geom_density(data = subset(flow_df, sample == "2i"), aes(x = Xist, ..scaled..), adjust = 0.8, 
               fill = "grey", alpha = 0.5, color = NA) +
  geom_density(data = subset(flow_df, sample != "2i"), aes(x = Xist, ..scaled.., color = sample), adjust = 0.8) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), limits = c(10, 200000),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Xist expression [a.u.]", y = "Cell density") +
  scale_color_manual(breaks = c("pLenti", "sgXist", "Library"), 
                     values = c("#12100B", "#92194C", "#CD2418"))

fix <- set_panel_size(facs_histogram, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("screen_density.pdf", fix, dpi = 300, useDingbats=FALSE) 


