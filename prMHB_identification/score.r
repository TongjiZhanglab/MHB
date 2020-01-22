# Ranking all potential allel-specific DHBs to identify potential regulatory MHBs
# input: identified stable maternal-specific MHBs and paternal-specific MHBs

mat_score <- read.csv("stable.maternalSpecific.DHBs.bed", header = FALSE, sep = "\t", stringsAsFactors = F)
pat_score <- read.csv("stable.paternalSpecific.DHBs.bed", header = FALSE, sep = "\t", stringsAsFactors = F)
colnames(mat_score) <- c("chrom", "start", "end", "DHB")
colnames(pat_score) <- c("chrom", "start", "end", "DHB")
mat_DHBs <- mat_score[,"DHB"]
pat_DHBs <- pat_score[,"DHB"]
rownames(mat_score) <- mat_DHBs
rownames(pat_score) <- pat_DHBs

asym_fold <- function(allele, df){
  colnames(df) <- c("chrom", "start", "end", "DHB", "oocyte", "sperm", "2cell.GG", "2cell.AG", "morula.GG", "morula.AG", "ICM.GG", "ICM.AG", "TE.GG", "TE.AG")  
  rownames(df) <- df[,"DHB"]
  DHBs <- rownames(df)
  if(allele=="M"){
    df_asym <- cbind( df[DHBs,"sperm"] / df[DHBs,"oocyte"], df[DHBs,"2cell.AG"] / df[DHBs,"2cell.GG"], df[DHBs,"morula.AG"] / df[DHBs,"morula.GG"], df[DHBs,"ICM.AG"] / df[DHBs,"ICM.GG"], df[DHBs,"TE.AG"] / df[DHBs,"TE.GG"])
  }else{
    df_asym <- cbind( df[DHBs,"oocyte"] / df[DHBs,"sperm"], df[DHBs,"2cell.GG"] / df[DHBs,"2cell.AG"], df[DHBs,"morula.GG"] / df[DHBs,"morula.AG"], df[DHBs,"ICM.GG"] / df[DHBs,"ICM.AG"], df[DHBs,"TE.GG"] / df[DHBs,"TE.AG"])
  }
  df_asym <- as.matrix(df_asym)
  rownames(df_asym) <- DHBs
  return(df_asym)
}

asym_score <- function(df, cutoff){
  df[df>=cutoff] <- 1007
  df[df<cutoff] <- 1
  df[df==1007] <- 0
  vec_sum <- as.vector(rowSums(df))
  return(vec_sum)
}

# ------------
#  H3K9me3
# ------------

mat_H3K9me3 <- read.csv("H3K9me3/stable.maternalSpecific.DHB_H3K9me3.txt", header = FALSE,sep = "\t", stringsAsFactors = F)
pat_H3K9me3 <- read.csv("H3K9me3/stable.paternalSpecific.DHB_H3K9me3.txt", header = FALSE,sep = "\t", stringsAsFactors = F)

mat_H3K9me3_fold <- asym_fold("M",mat_H3K9me3)
pat_H3K9me3_fold <- asym_fold("P",pat_H3K9me3)

mat_score$asym.H3K9me3 <- asym_score(mat_H3K9me3_fold[mat_DHBs,],0.5)
pat_score$asym.H3K9me3 <- asym_score(pat_H3K9me3_fold[pat_DHBs,],0.5)

# --------------------------
#  DNA methylation density
# --------------------------

mat_methyl <- read.csv("MethylationDensity/stable.maternalSpecific.DHB_MethylationDensity.txt", header = FALSE, sep = "\t", stringsAsFactors = F)
pat_methyl <- read.csv("MethylationDensity/stable.paternalSpecific.DHB_MethylationDensity.txt", header = FALSE, sep = "\t", stringsAsFactors = F)

mat_methyl_fold <- asym_fold("M",mat_methyl)
pat_methyl_fold <- asym_fold("P",pat_methyl)

mat_score$asym.methyl <- asym_score(mat_methyl_fold[mat_DHBs,],0.5)
pat_score$asym.methyl <- asym_score(pat_methyl_fold[pat_DHBs,],0.5)

# ------------------------
#   Regulatory elements
# -----------------------

# morula DHS (including allele-specific DHS)
mat_DHS <- read.csv("RegulatoryElements/DHS/stable.maternalSpecific.DHBs.morula.AG.DHS.bed", header = FALSE, sep = "\t", stringsAsFactors = F)
pat_DHS <- read.csv("RegulatoryElements/DHS/stable.paternalSpecific.DHBs.morula.GG.DHS.bed", header = FALSE, sep = "\t", stringsAsFactors = F)

mat_DHS_as <- unique(mat_DHS[which(mat_DHS[,5]>0 & mat_DHS[,6]>0),4])
mat_DHS_bi <- unique(mat_DHS[which(mat_DHS[,5]>0 & mat_DHS[,6]==0),4])
pat_DHS_as <- unique(pat_DHS[which(pat_DHS[,5]>0 & pat_DHS[,6]>0),4])
pat_DHS_bi <- unique(pat_DHS[which(pat_DHS[,5]>0 & pat_DHS[,6]==0),4])

mat_score$DHS <- 0
mat_score[mat_DHS_as,"DHS"] <- 2
mat_score[mat_DHS_bi,"DHS"] <- 1
pat_score$DHS <- 0
pat_score[pat_DHS_as,"DHS"] <- 2
pat_score[pat_DHS_bi,"DHS"] <- 1

# CTCF
mat_CTCF <- read.csv("RegulatoryElements/CTCF.lncRNA/stable.maternalSpecific.DHBs.CTCF.bed", header = F, sep = "\t", stringsAsFactors = F)
pat_CTCF <- read.csv("RegulatoryElements/CTCF.lncRNA/stable.paternalSpecific.DHBs.CTCF.bed", header = F, sep = "\t", stringsAsFactors = F)

mat_score$CTCF <- 0
mat_score[mat_CTCF[which(mat_CTCF[,5]>0),4],"CTCF"] <- 1
pat_score$CTCF <- 0
pat_score[pat_CTCF[which(pat_CTCF[,5]>0),4],"CTCF"] <- 1

# lncRNA
mat_lncRNA <- read.csv("RegulatoryElements/CTCF.lncRNA/stable.maternalSpecific.DHBs.lncRNA.bed", header = F, sep = "\t", stringsAsFactors = F)
pat_lncRNA <- read.csv("RegulatoryElements/CTCF.lncRNA/stable.paternalSpecific.DHBs.lncRNA.bed", header = F, sep = "\t", stringsAsFactors = F)

mat_score$lncRNA <- 0
mat_score[mat_lncRNA[which(mat_lncRNA[,5]>0),4],"lncRNA"] <- 1
pat_score$lncRNA <- 0
pat_score[pat_lncRNA[which(pat_lncRNA[,5]>0),4],"lncRNA"] <- 1

# methylated ZFP57 motif
mat_ZFP57 <- read.csv("RegulatoryElements/Methylated.ZFP57/stable.maternalSpecific.DHBs.ZFP57.morula.GG.methyl.bed", header = F, sep = "\t", stringsAsFactors = F)
pat_ZFP57 <- read.csv("RegulatoryElements/Methylated.ZFP57/stable.paternalSpecific.DHBs.ZFP57.morula.AG.methyl.bed", header = F, sep = "\t", stringsAsFactors = F)

mat_ZFP57_methylated <- unique(mat_ZFP57[which(mat_ZFP57[,5]>0.5),4])
pat_ZFP57_methylated <- unique(pat_ZFP57[which(pat_ZFP57[,5]>0.5),4])

mat_score$methylated.ZFP57 <- 0
mat_score[mat_ZFP57_methylated,"methylated.ZFP57"] <- 1
pat_score$methylated.ZFP57 <- 0
pat_score[pat_ZFP57_methylated,"methylated.ZFP57"] <- 1

# -------------
#  Expression
# -------------

# allele-specific gene expression in pre-implantation embryos (FC=2, stage = 2)
mat_expression <- read.csv("Expression/allelic.specific.2fold.2stage/stable.maternalSpecific.DHBs.350k.ASGene.bed", header = F,sep = "\t", stringsAsFactors = F)
pat_expression <- read.csv("Expression/allelic.specific.2fold.2stage/stable.paternalSpecific.DHBs.350k.ASGene.bed", header = F,sep = "\t", stringsAsFactors = F)
mat_score$asym.expression <- 0
mat_score[mat_expression[which(mat_expression[,5]==1),4],"asym.expression"] <- 1
pat_score$asym.expression <- 0
pat_score[pat_expression[which(pat_expression[,5]==1),4],"asym.expression"] <- 1

# allele-specific gene expression in pre-implantation embryos (FC=2, stage = 1)
mat_expression <- read.csv("Expression/allelic.specific.2fold.1stage/stable.maternalSpecific.DHBs.350k.ASGene.bed", header = F,sep = "\t", stringsAsFactors = F)
pat_expression <- read.csv("Expression/allelic.specific.2fold.1stage/stable.paternalSpecific.DHBs.350k.ASGene.bed", header = F,sep = "\t", stringsAsFactors = F)
mat_score$asym.expression <- 0
mat_score[mat_expression[which(mat_expression[,5]==1),4],"asym.expression"] <- 1
pat_score$asym.expression <- 0
pat_score[pat_expression[which(pat_expression[,5]==1),4],"asym.expression"] <- 1

# known imprinting genes
mat_KnownIG <- read.csv("Expression/Known.IG/stable.maternalSpecific.DHBs.350k.Known.IG.bed", header = F, sep = "\t", stringsAsFactors = F)
pat_KnownIG <- read.csv("Expression/Known.IG/stable.paternalSpecific.DHBs.350k.Known.IG.bed", header = F, sep = "\t", stringsAsFactors = F)
mat_score[mat_expression[which(mat_KnownIG[,5]==1),4],"asym.expression"] <- 1
pat_score[pat_expression[which(pat_KnownIG[,5]==1),4],"asym.expression"] <- 1

# ----------
#  ranking
# ---------

mat_score$asym.modification <- (mat_score$asym.H3K9me3 / 5 + mat_score$asym.methyl / 5) / 2
mat_score$regulatory.elements <- (mat_score$DHS + mat_score$CTCF + mat_score$lncRNA + mat_score$methylated.ZFP57) / 5
mat_score$score <- mat_score$asym.modification + mat_score$regulatory.elements + mat_score$asym.expression
# mat_score$score <- mat_score$asym.modification + mat_score$regulatory.elements

pat_score$asym.modification <- (pat_score$asym.H3K9me3 / 5 + pat_score$asym.methyl / 5) / 2
pat_score$regulatory.elements <- (pat_score$DHS + pat_score$CTCF + pat_score$lncRNA + pat_score$methylated.ZFP57) / 5
pat_score$score <- pat_score$asym.modification + pat_score$regulatory.elements + pat_score$asym.expression
# pat_score$score <- pat_score$asym.modification + pat_score$regulatory.elements

mat_score <- mat_score[order(mat_score$score,decreasing = TRUE),]
pat_score <- pat_score[order(pat_score$score,decreasing = TRUE),]

write.table(mat_score, 'stable.maternalSpecific.DHB.score.txt', row.names = F, col.names = T, sep = "\t", quote = F)
write.table(pat_score, 'stable.paternalSpecific.DHB.score.txt', row.names = F, col.names = T, sep = "\t", quote = F)

merge_score <- rbind(mat_score, pat_score)
merge_score <- merge_score[order(merge_score$score,decreasing = TRUE),]
write.table(merge_score, 'stable.alleleSpecific.DHB.score.txt', row.names = F, col.names = T, sep = "\t", quote = F)
