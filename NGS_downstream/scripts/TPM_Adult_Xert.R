#This script calculates TPM for Fig S3 (Xert expression in adult tissues)
library(Rsubread)
library(EnvStats)
library(Rgb)
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
soellner_dir <- args[1]
bauer_dir <- args[2]
wang_dir <- args[3]
gencode_xert <- args[4]
output_dir <- args[5]

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

setwd(soellner_dir)
soellner_temp = normalizePath(list.files(pattern="*.bam$", full.names = TRUE))

setwd(bauer_dir)
bauer_temp = normalizePath(list.files(pattern="*.bam$", full.names = TRUE))

setwd(wang_dir)
wang_temp = normalizePath(list.files(pattern="*.bam$", full.names = TRUE))

paired_temp <- c(bauer_temp, wang_temp)
  
feature_counts_soellner <- featureCounts(soellner_temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, 
                                         isPairedEnd = FALSE, GTF.featureType = "exon", strandSpecific = 0, 
                                         allowMultiOverlap = FALSE)

feature_counts_paired <- featureCounts(paired_temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, 
                                      isPairedEnd = TRUE, GTF.featureType = "exon", strandSpecific = 0, 
                                      allowMultiOverlap = FALSE)

counts_soellner <- data.frame(feature_counts_soellner$counts, feature_counts_soellner$annotation) %>%
  select(-Chr, -Start, -End, -Strand)

counts_paired <- data.frame(feature_counts_paired$counts, feature_counts_paired$annotation) %>%
  select(-Chr, -Start, -End, -Strand)


counts <- inner_join(counts_soellner, counts_paired)
names(counts) <- gsub(x = names(counts_paired), pattern = "\\.", replacement = "_")
names(counts) <- gsub(x = names(counts_paired), pattern = "_dedup_bam", replacement = "")


TPM <- counts %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length) %>% 
  inner_join(gene_names) %>% 
  select(-GeneID) %>%
  select(gene, everything())

setwd(output_dir)

write_delim(TPM, "S3E_adult_rnaseq_TPM.txt", delim = "\t")


TPM_long <- TPM %>%
  pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
  separate(sample, c("seq", "source", "tissue", "replicate"), sep = "_")


genes <- c("Xert")

tpm_plot <- TPM_long %>%
  filter(gene %in% genes) %>%
  ggplot(aes(x = tissue, y = tpm, color = factor(gene, levels = genes))) +
  facet_wrap(~factor(gene, levels = genes), scales = "free_y", ncol = 1) +
  geom_point(size = 0.6, alpha = 0.4) +
  stat_summary(fun = "mean", geom = "crossbar", lwd = 0.25, width = 0.5) +
  scale_color_manual(values = c("#00A69C")) +
  labs(y = "TPM") +
  ylim(0, 10)

pdf("S3E_Adult_Xert_TPM.pdf", useDingbats = FALSE, onefile = FALSE)
fix <- set_panel_size(tpm_plot, width = unit(4, "cm"), height = unit(1.25, "cm"))

print(grid.arrange(fix))

dev.off()


