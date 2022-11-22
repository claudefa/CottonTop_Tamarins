#!/bin/bash

#SBATCH -c 1
#SBATCH --mem-per-cpu 1G
#SBATCH --time=4:00:00
#SBATCH --array=1%10
#SBATCH --output=/projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/out/High.%A_%a.log
#SBATCH --job-name GL_High

dir="/projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${dir}Samples_5x | awk '{print $1}')
total=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | awk '{print $10}' | grep "/1" | wc -l) # total cout of derived alleles
totalhom=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | awk '{print $10}' | grep "1/1"| wc -l) # total cout of HOM derived alleles
totalLoadHigh=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep HIGH | awk '{print $10}' | grep "/1" | wc -l) # total load HIGH
RealizedLoadHigh=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep HIGH | awk '{print $10}' | grep "1/1" | wc -l ) # realized load HIGH
totalLoadMod=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep MODERATE | awk '{print $10}' | grep "/1" | wc -l) # total load HIGH
RealizedLoadMod=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep MODERATE | awk '{print $10}' | grep "1/1" | wc -l ) # realized load MODERATE
synonymous=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep synonymous_variant | awk '{print $10}' | grep "/1" | wc -l ) # total synonymous
synonymousHOM=$(vcftools --gzvcf ${dir}CTT_5x.ann.vcf.gz --indv $sample --stdout --recode --recode-INFO-all | grep synonymous_variant | awk '{print $10}' | grep "1/1" | wc -l ) # HOM synonymous
echo $sample | awk -v t=$total -v th=$totalhom -v tL=$totalLoadHigh -v RL=$RealizedLoadHigh -v tM=$totalLoadMod -v rM=$RealizedLoadMod -v S=$synonymous -v SH=$synonymousHOM '{print $1"\t"t"\t"th"\t"tL"\t"RL"\t"tM"\t"rM"\t"S"\t"SH}' > ${dir}Counts/${sample}_countsLoad.txt
