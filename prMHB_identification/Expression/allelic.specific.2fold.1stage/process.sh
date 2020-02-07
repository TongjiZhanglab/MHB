#!/bin/bash

# awk '{OFS="\t";FS="\t"}NR==FNR{a[$1]=$1}NR>FNR{if($5 in a)print $1,$2,$3,$4,$5,$6}' AS.FC2_1stages.genelist.GeneSymbol ~/data/Annotation/mm9/mm9.refseq.bed > AS.FC2_1stages.gene.genebody.bed
# awk '{OFS="\t";FS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' AS.FC2_1stages.gene.genebody.bed > AS.FC2_1stages.gene.tss.bed   


awk '{OFS="\t";FS="\t"}NR==FNR{a[$1]=$1}NR>FNR{if($4 in a){print $1,$2,$3,$4,$5}}' potential_IG_fc2_1stages_pre_RefSeq ~/data/Annotation/mm9/mm9.TSS.bed > AS.FC2_1stages.pre.genelist.tss
awk '{OFS="\t";FS="\t"}NR==FNR{a[$1]=$1}NR>FNR{if($4 in a){print $1,$2,$3,$4,$5}}' potential_IG_fc2_1stages_post_RefSeq ~/data/Annotation/mm9/mm9.TSS.bed > AS.FC2_1stages.post.genelist.tss 
