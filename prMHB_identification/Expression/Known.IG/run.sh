#!/bin/bash

awk '{OFS="\t"}NR==FNR{a[$1]=$0}NR>FNR{if($5 in a)print $1,$2,$3,$4,$5}' /mnt/Storage/home/yanghui/annotations/imprinted_genes_merged.euchr.GeneSymbol.txt ~/data/Annotation/mm9/mm9.TSS.bed > known.imprinting.gene.tss
sort -k1,1 -k2,2n known.imprinting.gene.tss | uniq > known.imprinting.gene.sort.tss
