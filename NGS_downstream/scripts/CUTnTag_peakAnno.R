#The script plots peak annotation for Supplementary Figure 2f from the d0 samples from the XXdXIC cell line
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(egg)
library(GenomicRanges)
library(ChIPseeker)
library(RColorBrewer)
library(extrafont)
loadfonts()

theme_set(theme_classic() +
  theme(legend.title = element_blank(), axis.title.y = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(),
                axis.text = element_text(size = 6, family = "Arial"), axis.title.x = element_text(size = 6, family = "Arial"),
                plot.title = element_blank(), 
                legend.text = element_text(size = 6, family = "Arial"), legend.key.size = unit(0.75,"line"),
                axis.ticks.y = element_blank()))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] #Supply path to BED files containing peaks (samples from XXdXIC/d0)
output_dir <- args[2]

setwd(input_dir)
samplefiles <- as.list(list.files(input_dir, pattern= "XXdXic.*d0_merged.bed", full.names=FALSE))
names(samplefiles) <- sapply(samplefiles, function(x) str_match(x, "(H[0-9,A-Z,a-z]*)")[,2])
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)

df_list <- list()

for (i in names(peakAnnoList)) {
  df <- data.frame(peakAnnoList[[i]]@anno) %>%
    mutate(annotation = ifelse(str_detect(annotation, "Intron"), "Intron",
                               ifelse(str_detect(annotation, "Downstream|Intergenic"), "Intergenic",
                               ifelse(str_detect(annotation, "Exon"), "Exon", annotation)))) %>%
    group_by(annotation) %>%
    summarize(n = n()) %>%
    mutate(percent = n/sum(n))
  df_list[[i]] <- df
}

total_peakAnno <- bind_rows(df_list, .id = "sample")

setwd(output_dir)

sample_order = c("H3K36me3", "H2AK119ub", "H3K27me3", "H3K9me3", "H3K4me1", "H3K27ac", "H3K4me3")
anno_order = c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic")

anno_plot <- total_peakAnno %>%
  ggplot(aes(y = factor(sample, levels = sample_order), x = percent, fill = factor(annotation, levels = anno_order))) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Fraction of Peaks") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_fill_manual(values = brewer.pal(7, "Set1")[c(1,2,3,4,5,7)])


pdf("Sup_Fig_2f_CUTnTag_peakAnno.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(anno_plot, width = unit(3, "cm"), height = unit(2, "cm"))

print(grid.arrange(fix))

dev.off()