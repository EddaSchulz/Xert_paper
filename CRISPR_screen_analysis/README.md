# CRISPR screen analysis

## Description
This folder contains all code that was used for alignment, counting and analysis of the CRISPR screen data in Gjaltema, Schwämmle 2021. 

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- MAGeCK (v0.5.9.3) collection of PYTHON scripts from "https://sourceforge.net/p/mageck/wiki/Home/"
- Python3 (v3.7.7) programming language from "https://www.python.org/"

R scripts can be run using R (v.3.6.3) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- extrafont (v0.17)
- gridExtra (v2.3)
- plyr (v1.8.6)
- tidyverse (v1.3.0)


## Reproduce analysis
Raw sequencing data from the CRISPRi screen can be retrieved from GSE167352. Data analysis and plotting can be performed using the master scripts in the "NGS_downstream/master/" directory. For all master scripts it is necessary to provide a path to "Xert_paper/NGS_alignment/" (with -p).

- (i)   "master_MAGeCK_run.sh": Conducts alignment and counting for the CRISPR screen data using MAGeCK count. Afterwards, it performs statistical analysis using MAGeCK mle and MAGeCK test. A path to a directory containing FASTQ files should be provided with -f. All other required files are stored in "/CRISPR_screen_analysis/files/".
- (ii)  "master_bootstrap_analysis.sh": Performs enrichment analysis in R using a bootstrapping approach and plots Figures 1G and S1M. An input count file is provided in "/CRISPR_screen_analysis/files/". Therefore, the script does not depend on "master_MAGeCK_run.sh".
- (iii) "master_mageck_analysis.sh": Uses the output of "master_MAGeCK_run.sh" to plot Figures 1D-F and S1K-L. All necessary input files are already provided in "CRISPR_screen_analysis/files/".
- (iv)  "master_screen_QC.sh": Plots correlation and guide distribution for Figures S1E and S1H-J. All necessary input files are already provided in "CRISPR_screen_analysis/files/". 


All scripts that are used by the master scripts are stored in "/CRISPR_screen_analysis/scripts/". Sample output for the master scripts is provided in "/CRISPR_screen_analysis/output/".
