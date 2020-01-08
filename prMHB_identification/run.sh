#!/bin/bash

path=`pwd`

# ---------------------------------------------
#  1. epigenetic modifications in stable 
#  allele-specific MHBs (DNA methylation amount and H3K9me3)
# ---------------------------------------------

# DNA methylation amount in stable allele-specific MHBs. Data source: raw data generated in the study.

cd ${path}/MethylationDensity

paste ${path}/stable.maternalSpecific.DHBs.bed \
stable.maternalSpecific.DHB_oocyte_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_sperm_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_2cell.GG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_2cell.AG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_morula.GG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_morula.AG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_ICM.GG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_ICM.AG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_TE.GG_MethylationDensity_siteprof1 \
stable.maternalSpecific.DHB_TE.AG_MethylationDensity_siteprof1 > stable.maternalSpecific.DHB_MethylationDensity.txt

paste ${path}/stable.paternalSpecific.DHBs.bed \
stable.paternalSpecific.DHB_oocyte_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_sperm_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_2cell.GG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_2cell.AG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_morula.GG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_morula.AG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_ICM.GG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_ICM.AG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_TE.GG_MethylationDensity_siteprof1 \
stable.paternalSpecific.DHB_TE.AG_MethylationDensity_siteprof1 > stable.paternalSpecific.DHB_MethylationDensity.txt

# H3K9me3 in stable allele-specific MHBs. Data source: raw data generated in the study.

cd ${path}/H3K9me3

paste ${path}/stable.maternalSpecific.DHBs.bed \
stable.maternalSpecific.DHB_oocyte_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_sperm_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_2cell.GG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_2cell.AG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_morula.GG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_morula.AG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_ICM.GG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_ICM.AG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_TE.GG_H3K9me3_siteprof1 \
stable.maternalSpecific.DHB_TE.AG_H3K9me3_siteprof1 > stable.maternalSpecific.DHB_H3K9me3.txt

paste ${path}/stable.paternalSpecific.DHBs.bed \
stable.paternalSpecific.DHB_oocyte_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_sperm_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_2cell.GG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_2cell.AG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_morula.GG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_morula.AG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_ICM.GG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_ICM.AG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_TE.GG_H3K9me3_siteprof1 \
stable.paternalSpecific.DHB_TE.AG_H3K9me3_siteprof1 > stable.paternalSpecific.DHB_H3K9me3.txt

# ---------------------------
#  2. Regulatory elements
# ---------------------------

# Overlapped with DHS. Data source: GSE92605
cd ${path}/RegulatoryElements/DHS
bedtools intersect -wao -a ${path}/stable.maternalSpecific.DHBs.bed -b morula.DHS.bed | \
bedtools intersect -wao -a - -b morula.AG.DHS.bed | cut -f 1-4,8,12 > stable.maternalSpecific.DHBs.morula.AG.DHS.bed
bedtools intersect -wao -a ${path}/stable.paternalSpecific.DHBs.bed -b morula.DHS.bed | \
bedtools intersect -wao -a - -b morula.GG.DHS.bed | cut -f 1-4,8,12 > stable.paternalSpecific.DHBs.morula.GG.DHS.bed

# Overlapped with CTCF or lncRNA. Data source: CTCF binding sites (GSE30206), lncRNA annotation (GENCODE)
cd ${path}/RegulatoryElements/CTCF.lncRNA
bedtools intersect -wao -a ${path}/stable.maternalSpecific.DHBs.bed -b CTCF.mm9.sorted.bed | cut -f 1-4,8 > stable.maternalSpecific.DHBs.CTCF.bed
bedtools intersect -wao -a ${path}/stable.paternalSpecific.DHBs.bed -b CTCF.mm9.sorted.bed | cut -f 1-4,8 > stable.paternalSpecific.DHBs.CTCF.bed
bedtools intersect -wao -a ${path}/stable.maternalSpecific.DHBs.bed -b gencode.mm9.lncRNA.sort.promoter | cut -f 1-4,9 > stable.maternalSpecific.DHBs.lncRNA.bed
bedtools intersect -wao -a ${path}/stable.paternalSpecific.DHBs.bed -b gencode.mm9.lncRNA.sort.promoter | cut -f 1-4,9 > stable.paternalSpecific.DHBs.lncRNA.bed

# overlapped with methylated ZFP57 motif TGCmCGC
cd ${path}/RegulatoryElements/Methylated.ZFP57 
python MotifLocate.py ${path}/stable.maternalSpecific.DHBs.bed stable.maternalSpecific.DHBs.MotifLocation.bed TGCCGC
python MotifLocate.py ${path}/stable.paternalSpecific.DHBs.bed stable.paternalSpecific.DHBs.MotifLocation.bed TGCCGC
cut -f 1-4 stable.maternalSpecific.DHBs.MotifLocation.bed > stable.maternalSpecific.DHBs.ZFP57.bed
cut -f 1-4 stable.paternalSpecific.DHBs.MotifLocation.bed > stable.paternalSpecific.DHBs.ZFP57.bed
bedtools intersect -wo -a stable.maternalSpecific.DHBs.ZFP57.bed -b /mnt/Storage/home/yuzhaowei/projects/imprinting/data/WGBS/morula.WT/morula.GG.1.sam.G.bed > stable.maternalSpecific.DHBs.ZFP57.morula.GG.methyl
bedtools intersect -wo -a stable.paternalSpecific.DHBs.ZFP57.bed -b /mnt/Storage/home/yuzhaowei/projects/imprinting/data/WGBS/morula.WT/morula.AG.1.sam.G.bed > stable.paternalSpecific.DHBs.ZFP57.morula.AG.methyl
python mcall_ave_methyl.py stable.maternalSpecific.DHBs.ZFP57.morula.GG.methyl stable.maternalSpecific.DHBs.ZFP57.morula.GG.methyl.bed
python mcall_ave_methyl.py stable.paternalSpecific.DHBs.ZFP57.morula.AG.methyl stable.paternalSpecific.DHBs.ZFP57.morula.AG.methyl.bed

# ------------------------------------
#  3. Allele-specific expressed genes
#     or known imprinting genes
# ------------------------------------

# allele-specific expressed genes within 350kb (fold change > 2 at least 2 stages)
cd ${path}/Expression/allelic.specific.2fold.2stage/
awk '{OFS=FS="\t"}{if($2>350000){print $1,$2-350000,$3+350000,$4}}' ${path}/stable.maternalSpecific.DHBs.bed > stable.maternalSpecific.DHBs.350k.bed
awk '{OFS=FS="\t"}{if($2>350000){print $1,$2-350000,$3+350000,$4}}' ${path}/stable.paternalSpecific.DHBs.bed > stable.paternalSpecific.DHBs.350k.bed
bedtools intersect -wao -a stable.maternalSpecific.DHBs.350k.bed -b AS.FC2_2stages.pre.genelist.tss | cut -f 1-4,9 | sort -k1,1 -k2,2n | uniq | awk '{OFS=FS="\t"}{if($5=="."){print $1,$2,$3,$4,0}else{print $1,$2,$3,$4,1}}' - | sort -k1,1 -k2,2n | uniq > stable.maternalSpecific.DHBs.350k.ASGene.bed
bedtools intersect -wao -a stable.paternalSpecific.DHBs.350k.bed -b AS.FC2_2stages.pre.genelist.tss | cut -f 1-4,9 | sort -k1,1 -k2,2n | uniq | awk '{OFS=FS="\t"}{if($5=="."){print $1,$2,$3,$4,0}else{print $1,$2,$3,$4,1}}' - | sort -k1,1 -k2,2n | uniq > stable.paternalSpecific.DHBs.350k.ASGene.bed

# known imprinting genes within 350kb
cd ${path}/Expression/Known.IG/
awk '{OFS=FS="\t"}{if($2>350000){print $1,$2-350000,$3+350000,$4}}' ${path}/stable.maternalSpecific.DHBs.bed > stable.maternalSpecific.DHBs.350k.bed
awk '{OFS=FS="\t"}{if($2>350000){print $1,$2-350000,$3+350000,$4}}' ${path}/stable.paternalSpecific.DHBs.bed > stable.paternalSpecific.DHBs.350k.bed
bedtools intersect -wao -a stable.maternalSpecific.DHBs.350k.bed -b known.imprinting.gene.sort.tss | cut -f 1-4,9 | sort -k1,1 -k2,2n | uniq | awk '{OFS=FS="\t"}{if($5=="."){print $1,$2,$3,$4,0}else{print $1,$2,$3,$4,1}}' - | sort -k1,1 -k2,2n | uniq > stable.maternalSpecific.DHBs.350k.Known.IG.bed
bedtools intersect -wao -a stable.paternalSpecific.DHBs.350k.bed -b known.imprinting.gene.sort.tss | cut -f 1-4,9 | sort -k1,1 -k2,2n | uniq | awk '{OFS=FS="\t"}{if($5=="."){print $1,$2,$3,$4,0}else{print $1,$2,$3,$4,1}}' - | sort -k1,1 -k2,2n | uniq > stable.paternalSpecific.DHBs.350k.Known.IG.bed

# ----------------------------
#  4. Calculating ranking score
# ----------------------------

Rscript score.r