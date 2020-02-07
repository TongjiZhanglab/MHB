#!/bin/bash

cp /mnt/Storage/home/yanghui/imprinting/data/RNA-Seq/t_data/potential_IG.fc2_2stage_pre_inhouseGA.RefSeq_GeneSymbol .
awk '{OFS="\t";FS="\t"}NR==FNR{a[$1]=$1}NR>FNR{if($4 in a){print $1,$2,$3,$4,$5}}' potential_IG.fc2_2stage_pre_inhouseGA.RefSeq_GeneSymbol ~/data/Annotation/mm9/mm9.TSS.bed > AS.FC2_2stages.pre.genelist.tss
