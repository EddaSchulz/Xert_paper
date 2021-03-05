# CRISPRi Screen Library Design

## Description
This folder contains all code that was used to design a composite set of candidate regulatory elements within the Xic for a CRISPR Screen. As an input the analysis uses ATAC-seq and STARR-seq data that was previously processed using the code in "/NGS_alignment/". 

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- MACS2 (v2.1.2) PYTHON script from "https://github.com/macs3-project/MACS"
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"

R scripts can be run using R (v.3.6.3) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- extrafont (v0.17)
- GenomicRanges (v1.38.0)
- gridExtra (v2.3)
- rtracklayer (v1.46.0)
- tidyverse (v1.3.0)
- UpSetR (v1.4.0)


## Reproduce analysis
The STARR-seq and ATAC-seq data (GSE167358) should be previously processed using the code within "/NGS_alignment/". Data analysis and plotting can be performed using the master scripts in the "/CRISPR_library_design/master/" directory. For all master scripts it is necessary to provide a path to "Xert_paper/NGS_alignment/".

- (i)   "master_library_design.sh":  Merges ATAC-seq and STARR-seq BAM files and calls peaks. Then combines all peaks within the Xic together with FANTOM5 data to form a composite peak set. For the final peak set used in the manuscript, REs with a length > 2kb were split up according to visual inspection of the ATAC-seq peaks. Guides were then retrieved using the Guidescan.com webtool. Requires paths containing ATAC-seq (with -a) and STARR-seq (with -s) BAM files.
- (ii)  "master_characterize_REs.sh": Characterizes candidate REs based on origin of peak, length and number of guides. Also plots this information. Input files are already provided under "/CRISPR_library_design/files/". Therefore, the script does not depend on "master_library_design.sh".


All scripts that are used by the master scripts are stored in "/CRISPR_library_design/scripts/". All files that are necessary to run the master scripts, other then the processed sequencing data, are stored in "/CRISPR_library_design/files/". Sample output of the master scripts is provided in "/CRISPR_library_design/output/". 
