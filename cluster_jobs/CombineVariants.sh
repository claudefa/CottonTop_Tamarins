#!/bin/bash

#Script to combine variants from GATK unified genotyper . 
#Ref genome
ref="/scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta"

#Directory
DIR="/scratch/devel/cfontser/CTT/VCF/"
OUTDIR=$DIR"/MergeVCFs/"
out=${OUTDIR}"out/";
qu=${OUTDIR}"qu/";

# create directories
mkdir -p $OUTDIR
mkdir -p $out
mkdir -p $qu

while read chrom;
do

	VCFs=""; while read line; do VCFs="$VCFs --variant ${DIR}/$chrom/${line}_${chrom}.vcf.gz"; done  < <(cat Samples);


	jobName=${qu}/CTT_${chrom}.combine.sh

	echo "java -Xmx4g -Djava.io.tmpdir=${TMPDIR} -jar /apps/GATK/3.5/GenomeAnalysisTK.jar -T CombineVariants -R ${ref} \
	$VCFs -o ${OUTDIR}/CTT_${chrom}.g.vcf.gz  -genotypeMergeOptions UNIQUIFY; tabix -p vcf ${OUTDIR}/CTT_${chrom}.g.vcf.gz " > ${jobName}

	chmod 755 $jobName
	echo $jobName

	python3 ~/submit.py  -i . -o ${out}/comb_${chrom}.out -e ${out}/comb_${chrom}.err -n comb_${chrom} -u 4 -t 1 -r lowprio -w 96:00:00 -c $jobName;
done < Chrom_autosomes
