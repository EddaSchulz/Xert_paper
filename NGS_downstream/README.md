# NGS downstream analysis

## Description
This folder contains all code that was used for analysis of sequencing data in Gjaltema, Schwämmle 2021. 

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- ChromHMM (v1.19) collection of JAVA scripts from "http://compbio.mit.edu/ChromHMM/"
- Deeptools2 (v3.4.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- diffReps (v1.55.6) PERL script from "https://github.com/shenlab-sinai/diffreps"
- Python3 (v3.7.7) programming language from "https://www.python.org/"

R scripts can be run using R (v.3.6.3) software ("https://cran.r-project.org/"). The following R libraries are required:
- ChIPseeker (v1.22.1)
- DESeq2 (v1.26.0)
- DiffBind (v2.6.6)
- egg (v0.4.5)
- EnvStats (v2.4.0)
- extrafont (v0.17)
- GenomicRanges (v1.38.0)
- gridExtra (v2.3)
- RColorBrewer (v1.1-2)
- Rgb (v1.6.1)
- Rsubread (v2.0.1)
- tidyverse (v1.3.0)
- TxDb.Mmusculus.UCSC.mm10.knownGene (v3.10.0)


## Reproduce analysis
The raw sequencing data (GSE167358) should be previously processed using the code within "/NGS_alignment/" (Output should not be renamed). Gene annotation (as GTF) has to be downloaded from https://www.gencodegenes.org/mouse/release_M25.html, appended with "/NGS_alignment/files/Xert.gtf" and stored as "/NGS_downstream/files/GENCODE_vM25_plus_Xert.gtf".
Data analysis and plotting can be performed using the master scripts in the "/NGS_downstream/master/" directory. For all master scripts it is necessary to provide a path to "Xert_paper/NGS_alignment/" (with -p) and to directories containing the various input files (that were previously generated with the code in "/NGS_alignment/".

- (i)   "master_ChromHMM.sh": Conducts ChromHMM analysis to retrieve putative RE states from CUT&Tag and ATAC-seq data. Returns colored BED files for use with the UCSC genome browser and plots the relative enrichment of different marks within different chromatin states. Requires paths to directories containing merged ATAC-seq (with -a) and CUT&Tag (with -c) BAM files (generated with "NGS_alignment/master/master_ATACseq.sh" and "NGS_alignment/master/master_CUTnTag.sh").
- (ii)   "master_CTCFtrack.sh: Generates a BED track containing directional CTCF binding sites (CBSs). As an input the script uses a TSV file that was generated from the FIMO webtool (https://meme-suite.org/meme/doc/fimo.html) using peaks generated from a male ESC line (GSE30203) with "/NGS_alignment/master/unpaired_ChIPseq.sh" and a CTCF motif file from JASPAR (http://jaspar.genereg.net/matrix/MA0139.1/). The TSV file is provided under "/NGS_downstream/files/fimo_ctcf.tsv". 
- (iii)  "master_CUTnTag_QC.sh": Performs several quality control analysis on the CUT&Tag data: It calculates pearson correlation between CUT&Tag samples, compares CUT&Tag data to native ChIP-seq in the parental TX1072 cell line via PCA and annotates peaks to gene body regions using ChIPseeker. Requires paths containing replicate CUT&Tag BAM files (with -b), merged CUT&Tag BAM files (with -m), merged CUT&Tag peak files (with -n) and replicate BAM files from Zylicz et al. 2019 ChIP-seq data (with -z). The input data should be generated previously with "NGS_alignment/master/master_CUTnTag.sh" and "NGS_alignment/master/master_Zylicz_ChIPseq.sh".
- (iv)   "master_diffEnrichment.sh": Performs DiffBind for ATAC-seq, H3K4me1, H3K4me3 and H3K27ac, as well as diffReps for H3K9me3 in various comparisons (X-dosage and differentiation-wise). Furthermore, it creates colored BED-files for use in UCSC. Requires ATAC-seq (-a for peaks/ -t for BAM files) and CUT&Tag (-c for peaks/ -b for BAM files) data that was previously processed using "/NGS_alignment/master_ATACseq.sh" and "/NGS_alignment/master_CUTnTag.sh".
- (v)   "master_karyotyping.sh": Compares chromosomewide reads from all novel cell lines in Gjaltema, Schwämmle 2021 to a TX1072 XX wildtype control in order to identify deletions or duplications. An input count table is already provided at "/NGS_downstream/files/karyotype_counts.txt" (the reads were previously quantified using multiBamSummary on the usegalaxy.eu platform). 
- (vi) "master_TTseq_downstream.sh": Counts reads within RNA-seq and TT-seq data. Then calculates TPM and identifies differentially expressed genes between the XXdXic and XO cell lines. Lastly, plots expression data for lncRNA genes in the Xic as a heatmap and a lineplot. Requires paths containing RNA-seq (with -r) and TT-seq (with -t) BAM files (generated with "NGS_alignment/master/master_RNAseq.sh" and "NGS_alignment/master/master_TTseq.sh"). 
- (vii)  "master_Zhang_RNAseq_downstream.sh": Counts reads and calculates TPM from published embryo RNA-seq data. Afterwards, it plots expression of several genes (Xert, Jpx, Ftx, Rnf12, Oct4 and Otx2). Requires path to directory containing processed embryo RNA-seq data (with -b). The input data should be previously generated with "NGS_alignment/master/master_Zhang_RNAseq.sh" 


All scripts that are used by the master scripts are stored in "/NGS_downstream/scripts/". All files that are necessary to run the master scripts, other then the GENCODE annotation and the processed sequencing data, are stored in "/NGS_downstream/files/". Sample output of the master scripts is provided in "/NGS_downstream/output/". 
