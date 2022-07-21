#!/bin/bash

# Set directories
DIR="/scratch/devel/cfontser/CTT/"
IN=${DIR}"BAMs/"
OUTDIR=${DIR}"VCF/"

mkdir -p $OUTDIR
# Set output
qu=${OUTDIR}"qu"
out=${OUTDIR}"out"

# Create directories

mkdir -p $qu
mkdir -p $out

# Reference

#REF
ASSEMBLY="/scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta"

module purge;module load gcc/6.3.0 popt/1.16 nlopt/2.5.0 R/3.5.0 snpad PYTHON/3.6.3 tabix

while read chrom;
do
while read line;
do
        sample=`echo $line | awk '{print $1}'`

        bam=${IN}${sample}.SagMidas.realigned.bam
	snpad=${OUTDIR}/${chrom}/${sample}_${chrom}.snpAD
	vcf=${OUTDIR}/${chrom}/${sample}_${chrom}.vcf
	mkdir -p ${OUTDIR}/${chrom}/
        jobName=${qu}/2VCF.${sample}_${chrom}.sh

	echo "Bam2snpAD -r $chrom -f $ASSEMBLY -Q 30 $bam  > $snpad; \
		snpAD -c 12 -o ${OUTDIR}/${chrom}/${sample}.priors.txt -O ${OUTDIR}/${chrom}/${sample}.errors.txt $snpad ; \
		snpADCall -N 1 -e ${OUTDIR}/${chrom}/${sample}.errors.txt -p \"\`cat ${OUTDIR}/${chrom}/${sample}.priors.txt\`\" $snpad > $vcf; \
		bgzip $vcf; tabix -p vcf ${vcf}.gz; rm $snpad; rm ${OUTDIR}/${chrom}/${sample}.priors.txt; rm ${OUTDIR}/${chrom}/${sample}.errors.txt" > ${jobName}
        chmod 755 $jobName
        echo $jobName
	python3 ~/submit.py  -c $jobName -n vcf${sample} -o ${out}/vcf${sample}_${chrom}.out -e  ${out}/vcf${sample}_${chrom}.err -u 8 -r lowprio -w "28:59:00"

done < Samples

done < Chrom_test
