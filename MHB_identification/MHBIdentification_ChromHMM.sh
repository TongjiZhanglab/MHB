#!/bin/bash

# ----------------------------------------------------------------------
#  MHB identification based on DNA methylation amount and 
#  H3K9me3 signals in mouse pre-implantation embryos
#  Input: big wiggle file for DNA methylation amount and H3K9me3 signals 
# ---------------------------------------------------------------------

dirPATH=`pwd`

chromHMM_Overall(){
	for sample in 2cell.GG 2cell.AG morula.GG morula.AG ICM.GG ICM.AG TE.GG TE.AG
	do

		if [ ! -d ${dirPATH}/chromHMM/${sample} ];then mkdir -p ${dirPATH}/chromHMM/${sample};fi;cd ${dirPATH}/chromHMM/${sample}
		# Step1. Construct config file for generating ChromHMM's binarize signal
		echo -e "${sample}\tmethylDensity\tmm9.1kbB_10bpS.methylDensity.${sample}.wt.bw" > config_bw.txt	# DNA methylation amount, data source: GSE97778
		echo -e "${sample}\tH3K9me3\t${sample}.wt.H3K9me3.merged.32M.macs2.bw" >> config_bw.txt				# H3K9me3 signals, data source: GSE97778

		# Step2. Rebin and combine from signal files
		python ${dirPATH}/chromhmm_signalize.py rebin --binsize 200 --binned_dir binned_bw_200 config_bw.txt mm9
		rm ${dirPATH}/chromHMM/${sample}/binned_bw_200/*-chrM.binned ${dirPATH}/chromHMM/${sample}/binned_bw_200/*-chrX.binned ${dirPATH}/chromHMM/${sample}/binned_bw_200/*-chrY.binned
		python ${dirPATH}/chromhmm_signalize.py combine --binsize 200 --binned_dir binned_bw_200 --combined_dir signal_bw_200 config_bw.txt

		# step3. Calculate balance parameter (lambda) for each chromosome
		if [ ! -d ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_lambda ];then mkdir -p ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_lambda;fi;cd ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_lambda
		for num in $(seq 1 19)
		do
			Rscript ${dirPATH}/balancedLambdaForEachChromosome.r chr${num} ${dirPATH}/chromHMM/${sample}/signal_bw_200/${sample}_chr${num}_200_signal.txt ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_lambda/lambdaForEachChromosome.txt
		done

		# Step4. Calculated average lambda of chr1-19.
		lambda_overall=`cat ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_lambda/lambdaForEachChromosome.txt | awk 'BEGIN{FS=OFS="\t";T=0}{T+=$2}END{print T/NR}'`

		# Step5. Balance DNA methylation amount and H3K9me3 according to overall lambda
		if [ ! -d ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_OverallLambda ];then mkdir -p ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_OverallLambda;fi
		for num in $(seq 1 19)
		do
			Rscript ${dirPATH}/balanceForEachChromosomeAccordingToLambda.r ${lambda_overall} ${dirPATH}/chromHMM/${sample}/signal_bw_200/${sample}_chr${num}_200_signal.txt ${dirPATH}/chromHMM/${sample}/signal_bw_200_balanced_OverallLambda/${sample}_chr${num}_200_signal.txt
		done

		# Step6. MHBs identification
		for foldthreshold in 0
		do
			for pvalue in 0.000005
			do
				cd ${dirPATH}/chromHMM/${sample}
				java -Xmx50000m -jar ChromHMM.jar BinarizeSignal -f ${foldthreshold} -p ${pvalue} signal_bw_200_balanced_OverallLambda binarized_bw_200_balanced_OverallLambda_f${foldthreshold}_p${pvalue}
				java -Xmx50000m -jar ChromHMM.jar LearnModel -b 200 -p 24 binarized_bw_200_balanced_OverallLambda_f${foldthreshold}_p${pvalue} output_bw_200_balanced_OverallLambda_f${foldthreshold}_p${pvalue} 4 mm9
				cd ${dirPATH}/chromHMM/${sample}/output_bw_200_balanced_OverallLambda_f${foldthreshold}_p${pvalue}
				State=$(cat emissions_4.txt | awk 'BEGIN{FS=OFS="\t";S="E1";MAX=-1;}{if($2+$3>MAX){MAX=$2+$3;S=$1}}END{print "E"S}')
				grep -w ${State} ${sample}_4_segments.bed | sort -k1,1 -k2,2n | mergeBed -i - -d 3000 | awk 'BEGIN{FS=OFS="\t"}{if($3-$2>=600){print $1,$2,$3}}' > balanced_OverallLambda_f${foldthreshold}_p${pvalue}_dist3k_MHB_len600_${sample}.bed
			done
		done
	done
}
