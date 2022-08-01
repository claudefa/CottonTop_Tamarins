#!/bin/bash

bed=$1;
sample=$2;
if [ -z "$bed" ] || [ -z "$sample" ];then printf "### Calculate Het by windows ###\n\nUSAGE:\n\nbash calculateHet_windows.sh 100000 Sample \n\n";exit;fi
mkdir -p Het_windows_real/${sample}/

while read line;
do 
	region=$(echo $line | awk '{print $1":"$2"-"$3}');
	read chrom start end <<<$(echo $line);
	vcf="/scratch/devel/cfontser/CTT/VCF/MergeVCFs/CTT_${chrom}.g.vcf.gz"
	callable=$(tabix -h $vcf $region | vcftools --vcf - --indv $sample --max-missing 1 --stdout --recode --recode-INFO-all --minDP 3 --maxDP 30 --minQ 30 | grep -v '#' | wc -l ) 
	if [ $callable -gt 0 ];
	then
		tabix -h $vcf $region | vcftools --vcf - --indv $sample --max-missing 1 --stdout --recode --recode-INFO-all --minDP 3 --maxDP 30 --minQ 30 | vcfhetcount | tail -n 1 | awk -v l="$callable" -v chrom="$chrom" -v start="$start" '{print chrom"\t"start+50000"\t"$1"\t"l"\t"$1/l}'
	else
		echo $line | awk -v l=0 -v het=NA '{print $1"\t"$2+50000"\t"l"\t"l"\t"het}'
	fi
	callable=''
done <  Windows/${bed}.bed > Het_windows_real/${sample}/${sample}_${bed##*/}_wind.het
