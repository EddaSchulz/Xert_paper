# Generates a BED file with directional CTCF binding sites (data from Stadler 2011)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] 
output_dir <- args[2]

setwd(input_dir)

fimo_tsv <- read.delim("fimo_ctcf.tsv")

setwd(output_dir)

fimo <- fimo_tsv %>%
  separate(sequence_name, c("chr", "seq"), sep = ":") %>%
  separate(seq, c("seq_start", "seq_end"), sep = "-") %>%
  mutate(start = as.numeric(seq_start) + as.numeric(start) - 1, end = as.numeric(seq_start) + as.numeric(stop) - 1) %>%
  transmute(chr = chr, start = start, end = end, name = ".", score = 1000, strand = strand, thickStart = start, thickEnd = end,
            rgb = ifelse(strand == "+", "255,000,000", "000,000,255")) %>%
  arrange(chr, start)


write_delim(fimo, "S7A_FIMO_CBS.bed", delim = "\t", col_names = FALSE)
