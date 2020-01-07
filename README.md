# Source codes for MHB study



Allele-specific DNA methylation and H3K9me3 blocks are critical for early embryo development



Our study comprehensively analyzed the association between epigenetic features and DNA methylation maintenance in mouse pre-implantation embryos, and revealed that H3K9me3 can facilitate DNA methylatin maintenance in genome loci with the co-existence of DNA methylation and H3K9me3, termed as DNA methylation and H3K9me3 blocks (MHBs). We further profiled parental MHBs in pre-implantation embryos separately and identified 502 stable allele-specific MHBs, including 17 out of 20 known gamete imprinting control regions. We experimentally validated that H3K9me3 is essential for the repressive function of such allele-specific MHBs. Moreover, we identified 31 potential regulatory MHBs may regulate gene transcription like imprinting control regions, and two of them are validated to be critical for embryonic development. Taken together, our study revealed the widespread existence of allele-specific MHBs and their transcription regulation functions, largely extended the scope of allele-specific regulation in mammalian pre-implantation embryos. The study has been warpped up to submit.



Here, we provided main source codes for the study.

#### 1. Identification of DNA methylation and H3K9me3 blocks

---------------------------------------------------------------------------------------------------------------------------------------------------------------

> To identify MHBs in pre-implantation embryos, we first divided the mouse genome into 200-bp consecutive bins. We then calculated the average H3K9me3 signals and DNA methylation amount on such regions. Given that the significant difference of values of H3K9me3 signals and DNA methylation amount, we transformed H3K9me3 signals by multiplying a parameter which was defined as the fold of the mean of DNA methylation amount and the mean of H3K9me3 signals. We then applied ChromHMM on DNA methylation amount and transformed H3K9me3 signals at such bins and identified bins with high DNA methylation and H3K9me3 simultaneously as MHBs. Neighbor MHBs no more than 3000 bp apart were merged and MHBs shorter than 600 bp were removed. 

+ The source code for identification of DNA methylatin and H3K9me3 blocks is available at 	[MHB_identification]().

#### 2. Identificaiton of potential regulatory allele-specific MHBs
------------------------------------------------------------------------------------------------------------------------------

> To identify potential regulatory allele-specific MHBs with similar transcription regulation function to known ICR, we calculated a score for all identified allele-specific MHBs based on some features of known ICRs. We scored allele-specific MHBs according to three categories of features, epigenetic modification imbalance, regulatory elements and surrounding allele-specific expressed genes. For epigenetic modification imbalance, we evaluated the level of imbalance of DNA methylation amount and H3K9me3 signals at gametes, 2-cell, morula, ICM and TE stages for each MHB. The score of epigenetic modification for each MHB was denoted by the number of stages where the given MHB showed imbalance of the epigenetic modification between GG and AG embryos (for maternal specific MHB, it meant GG / AG > 2 and for paternal specific MHB, it meant AG / GG > 2). For regulatory elements, we scored MHBs using the following four subcategories: (1). MHBs overlapped with allele-specific DNase I hypersensitive sites (DHSs) or bi-allelic DHSs were scored to 2 or 1 points, respectively. DNase-seq data of GG and AG embryos were from Inoue et al.’s study; (2) MHBs overlapped with CTCF binding sites in mESC were scored to 1 point, ChIP-seq data of CTCF in mESC was from Stadler et al.’s study; (3) MHBs overlapped with lncRNA TSSs were scored to 1 point, the annotation file of lncRNA was from GENCODE; (4) MHBs overlapped with methylated Zfp57 motif (TGCmCGC) (DNA methylation level on the CpG of hexanucleotide > 0.5) were scored to 1 point. For surrounding allele-specific expressed genes, we first defined genes with expression imbalance between GG and AG embryos (fold change of FPKM values > 2) in at least two pre-implantation stages as allele-specific expressed genes. MHBs with the existence of allele-specific expressed genes (fold change of FPKM values > 2) or known imprinting genes within 350kb were score to 1 point. We normalized three categories of features to be the same weight in rank score, i.e. the total points of three categories of features was set to 1, respectively. Then, we merged points of all three categories of features to obtain a rank score for all MHBs. Scores of all identified allele-specific MHBs were annotated in Supplementary Table 3 in the paper. MHBs had equal or higher score than known ICR were identified as potential regulatory allele-specific MHBs.

+ The source code for identification of potential regulatory allele-specific MHBs is available at [prMHB_identification]().

  
