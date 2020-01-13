# Source codes for MHB study



Allele-specific DNA methylation and H3K9me3 blocks are critical for early embryo development



Our study profiled allele-specific MHBs in mouse pre-implantation GG / AG embryos. Among more than 500 stable allele-specific MHBs, we identified dozens of potential regulatory MHBs with similar features to known imprinting control regions that may regulate gene transcription and experimentally validated two of them. The study has been warpped up to submit.



Here, we provided main source codes for the study.

#### 1. Identification of DNA methylation and H3K9me3 blocks

---------------------------------------------------------------------------------------------------------------------------------------------------------------

> To characterize the features of MHBs in pre-implantation embryos, we identified 6,002, 6,197, 6,622 and 8,108 MHBs at 2-cell, 8-cell, morula and ICM stages separately, by applying ChromHMM on DNA methylation amount and H3K9me3 signals.

+ The source code for identification of MHBs is available at	[MHB_identification]().

#### 2. Identificaiton of potential regulatory MHBs
------------------------------------------------------------------------------------------------------------------------------

> We calculated a score for each stable allele-specific MHB based on the features of known ICRs, including the allele-imbalance of DNA methylation amount and H3K9me3 signals, the overlapping with DNase I hypersensitive sites, Zpf57 binding motifs, CTCF binding sites and lncRNA genes. 

+ The source code for identification of potential regulatory allele-specific MHBs is available at [prMHB_identification]().

  
