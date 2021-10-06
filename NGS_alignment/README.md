# NGS Read Alignment and Data Processing

## Description
This folder contains the code for alignment and data processing of the following sequencing data:

From the present study: 

ATAC-seq: ATAC-seq was performed in the TX1072 XXΔXic (XX line carrying a heterozygous deletion of the Xic: chrX:103182701-103955531 on the B6 allele) and TX1072 XO (X chromosome from Cast background) cell lines. The data was collected in replicates from cells in 2i+LIF medium (day 0), as well at day 2 and day 4 of differentiation (via 2i+LIF-withdrawal).

CUT&Tag: CUT&Tag targeting H3K4me3, H3K27ac, H3K4me1, H3K9me3, H3K27me3, H3K36me3 and H2AK119ub was performed in the TX1072 XXΔXic (XX line carrying a heterozygous deletion of the Xic: chrX:103182701-103955531 on the B6 allele) and TX1072 XO (X chromosome from Cast background) cell lines. The data was collected in replicates from cells in 2i+LIF medium (day 0), as well at day 2 and day 4 of differentiation (via 2i+LIF-withdrawal).

pA-RNA-seq: RNA-seq enriched for polyadenylated fragments was performed in differentiating TX1072 XX mESCs (2 days of 2i+LIF-withdrawal) using the TruSeq RNA Sample Preparation Kit v2 (Illumina).

RNA-seq: RNA-seq was performed in the TX1072 XXΔXic (XX line carrying a heterozygous deletion of the Xic: chrX:103182701-103955531 on the B6 allele) and TX1072 XO (X chromosome from Cast background) cell lines as an input for the TT-seq protocol. The data was collected in replicates from cells in 2i+LIF medium (day 0), as well at day 2 and day 4 of differentiation (via 2i+LIF-withdrawal).

STARR-seq: STARR-seq was performed in the TX1072 XX and TX1072 XO cell lines using the BAC clones of the RPCI-23 Female (C57BL/6J) Mouse BAC Library (RP23-106C4, RP23-11P22, RP23-423B1, RP23-273N4, RP23-71K8) as the input DNA source. The data was collected in 3 replicates from cells in 2i+LIF medium (day 0), as well at day 2 of differentiation (via 2i+LIF-withdrawal).

TT-seq: TT-seq was performed in the TX1072 XXΔXic (XX line carrying a heterozygous deletion of the Xic: chrX:103182701-103955531 on the B6 allele) and TX1072 XO (X chromosome from Cast background) cell lines. The data was collected in replicates from cells in 2i+LIF medium (day 0), as well at day 2 and day 4 of differentiation (via 2i+LIF-withdrawal).

From published data:

ATAC-seq (Bleckwehl et al., 2021): Paired-End ATAC-seq data in XY mESCs and epiblast-like cells was retrieved from GSE155089.

ATAC-seq (Guo et al., 2017): Paired-End ATAC-seq data in XX hESCs was retrieved from GSE92280.

ChIP-seq (Bleckwehl et al., 2021): Single-End ChIP-seq data for H3K27ac and H3K9me3 in XY mESCs and epiblast-like cells was retrieved from GSE155089.

ChIP-seq (Buecker et al., 2014): Single-End ChIP-seq data for OCT4 and OTX2 in XX mESCs and epiblast-like cells, as well as the corresponding input samples, was retrieved from GSE56098.

ChIP-seq (Chronis et al., 2017): Single-End ChIP-seq data for C-MYC, ESRRB, KLF4, NANOG, OCT4 and SOX2 in XY mESCs was retrieved from GSE90893.

ChIP-seq (Gontan et al., 2012): Single-End ChIP-seq data for REX1 in XX mESCs was retrieved from GSE36417.

ChIP-seq (Li et al., 2019): Single-End ChIP-seq data for SMAD2/3 in XY hESCs, as well as the corresponding input samples, was retrieved from GSE109524.

ChIP-seq (Rada-Iglesias et al., 2011): Single-End ChIP-seq data for H3K4me1, H3K4me3, H3K27ac and P300 in XX hESCs, as well as the corresponding input samples, was retrieved from GSE28874.

ChIP-seq (Stadler et al., 2011): Single-End ChIP-seq data for CTCF in mESCs, as well as the corresponding input samples, was retrieved from GSE70486.

ChIP-seq (Tu et al. 2016): Single-End ChIP-seq data for PRDM14 in mESCs was retrieved from GSE71675.

ChIP-seq (Wang et al., 2017): Paired-End ChIP-seq data for SMAD2/3 and TCF3 in mESCs and embryoid bodies, as well as the corresponding input samples, was retrieved from GSE30203.

ChIP-seq (Zyclicz et al. 2019): Paired-End ChIP-seq data for H3K4me1, H3K4me3, H3K27ac, H3K27me3 and H2AK119ub in TX1072 XX mESCs was retrieved from GSE116990.

GRO-seq (Wang et al., 2015): Single-End GRO-seq data in XY ESCs was retrieved from GSE54471.

RNA-seq (Bauer et al., 2021): Paired-End RNA-seq data in XX NPCs was retrieved from GSE157448.

RNA-seq (Söllner et al., 2017): Single-End RNA-seq data in sex-mixed adult mouse tissues was retrieved from E-MTAB-6081.

RNA-seq (Wang et al., 2019): Paired-End RNA-seq data in XX MEFs was retrieved from GSE116413.

RNA-seq (Zhang et al., 2018):  Paired-End RNA-seq data from E3.5, E4.0, E5.5, E6.5 and E7.5 mouse embryos was retrieved from GSE76505.

scRNA-seq (Deng et al., 2014): Single-End scRNA-seq data from sex-mixed early mouse embryos (Zygote-16-cell stages) was retrieved from GSE45719.

scRNA-seq (Petropoulos et al., 2016): Single-End scRNA-seq data from sex-mixed human embryos (E3-E7 stages) was retrieved from E-MTAB-3929.

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- Bowtie2 (v2.3.5.1) aligner from "https://github.com/BenLangmead/bowtie2"
- Bowtie (v1.2.2) aligner from "http://bowtie-bio.sourceforge.net/index.shtml"
- Deeptools2 (v3.4.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- MACS2 (v2.1.2) PYTHON script from "https://github.com/macs3-project/MACS"
- Picard tools (v2.18.25) collection of JAVA scripts from "https://broadinstitute.github.io/picard/"
- Python3 (v3.7.7) programming language from "https://www.python.org/"
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- STAR (v2.5.4a) aligner from "https://github.com/alexdobin/STAR"
- TrimGalore (v0.6.4) PERL script from "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"


## Reproduce analysis
The alignment and data processing for the respective NGS data types can be performed using the master scripts stored in the "NGS_alignment/master/" directory. For all master scripts the respective FASTQ files have to be downloaded from GEO and renamed according to "NGS_alignment/files/rename_GSM.txt". The mm10 reference genome needs to be downloaded from ncbi and stored as "/NGS_alignment/files/mm10.fa" (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/). The B6/Cast-masked mm10 genome has to be prepared according to https://github.com/EddaSchulz/Pacini_paper and stored as "/NGS_alignment/files/N_masked_B6_Cast.fa". Gene annotation (as GTF) has to be downloaded from https://www.gencodegenes.org/mouse/release_M25.html, appended with "/NGS_alignment/files/Xert.gtf" and stored as "/NGS_alignment/files/GENCODE_vM25_plus_Xert.gtf". For human data, the hg38 reference genome needs to be downloaded from ncbi and stored as "/NGS_alignment/files/hg38.fa" (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)

For all scripts, it is necessary to provide a path to "Xert_paper/NGS_alignment/" (with -p) and the directory containing the FASTQ files (with -f). A work directory can be supplied with -d.

- (i)   "master_ATACseq.sh": Builds masked mm10 genome, aligns and processes ATAC-seq data, calls peaks and creates coverage tracks.
- (ii)  "master_CUTnTag.sh": Builds masked mm10 genome, aligns and processes CUT&Tag data, calls peaks and creates coverage tracks.
- (iii) "master_pA-RNAseq.sh": Builds masked mm10 genome, aligns and processes pA-RNAseq data and creates coverage tracks.
- (iv)  "master_RNAseq.sh": Builds masked mm10 genome, aligns and processes RNA-seq data.
- (v) "master_STARRseq.sh": Builds mm10 genome, aligns and processes STARR-seq data.
- (vi) "master_TTseq.sh": Builds masked mm10 genome, aligns and processes TT-seq data and creates coverage tracks.
- (vii) "master_unpaired_ChIPseq.sh": Builds mm10 genome, aligns and processes ChIP-seq data from Buecker et al. 2014 and Stadler et al. 2011, calls peaks and creates coverage tracks.
- (viii) "master_Wang_ChIPseq.sh": Builds mm10 genome, aligns and processes ChIPseq data from Wang et al. 2017, Calls peaks and creates coverage tracks.
- (ix) "master_Zhang_RNAseq.sh": Builds mm10 genome, aligns and processes embryo RNA-seq data from Zhang et al. 2018.
- (x) "master_Zylicz_ChIPseq.sh": Builds masked mm10 genome, aligns and processes ChIPseq data from Zylicz et al. 2019.

All scripts that are used by the master scripts are stored in "/NGS_alignment/scripts/". All files that are necessary to run the master scripts, other than those mentioned above (FASTQ files, genome annotation files), are stored in "/NGS_alignment/files/".
