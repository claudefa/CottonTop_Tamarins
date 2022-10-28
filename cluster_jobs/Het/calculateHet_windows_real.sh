#!/bin/bash

bed=$1;
sample=$2;
id=$3
if [ -z "$bed" ] || [ -z "$sample" ] || [ -z "$id" ] ;then printf "### Calculate Het by windows ###\n\nUSAGE:\n\nbash calculateHet_windows.sh 100000 Sample ID_vcf\n\n";exit;fi
mkdir -p /projects/mjolnir1/people/qvw641/CottonTop/VCF/Het/Het_windows_real/${sample}/

while read line;
do 
	region=$(echo $line | awk '{print $1":"$2"-"$3}');
	read chrom start end <<<$(echo $line);
	vcf="/projects/mjolnir1/people/qvw641/CottonTop/VCF/MergeVCFs/CTT_${chrom}.g.vcf.gz"
	callable=$(tabix -h $vcf $region | vcftools --vcf - --indv $id --max-missing 1 --stdout --recode --recode-INFO-all --minDP 3 --maxDP 30 --minQ 30 | grep -v '#' | wc -l ) 
	if [ $callable -gt 0 ];
	then
		tabix -h $vcf $region | vcftools --vcf - --indv $id --max-missing 1 --stdout --recode --recode-INFO-all --minDP 3 --maxDP 30 --minQ 30 | grep 0/1 | wc -l  | awk -v l="$callable" -v chrom="$chrom" -v start="$start" '{print chrom"\t"start+50000"\t"$1"\t"l"\t"$1/l}'
	else
		echo $line | awk -v l=0 -v het=NA '{print $1"\t"$2+50000"\t"l"\t"l"\t"het}'
	fi
	callable=''
done <  /home/qvw641/CottonTop_Tamarins/cluster_jobs/Het/Windows/${bed}.bed > /projects/mjolnir1/people/qvw641/CottonTop/VCF/Het/Het_windows_real/${sample}/${sample}_${bed##*/}_wind.het
