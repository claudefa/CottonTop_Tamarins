#!/bin/bash

# repeat analysis

# Set directories
DIR="/projects/mjolnir1/people/qvw641/CottonTop/"
IN=${DIR}/"BAMs/"
OUTDIR=${DIR}"VCF/"

mkdir -p $OUTDIR
# Set output
qu=${OUTDIR}"qu"
out=${OUTDIR}"out"

# Create directories

mkdir -p $qu
mkdir -p $out

# Reference
ASSEMBLY="/projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta"

module load snpAD/0.3.10

while read line
do
	chrom=`echo $line | awk '{print $2}'`
	sample=`echo $line | awk '{print $1}'`
	jobName=$qu/TXT_${chrom}_$sample.sh

	bam=${IN}${sample}.SagMidas.autosome.bam
	snpad=${OUTDIR}/${chrom}/${sample}_${chrom}.snpAD
	mkdir -p ${OUTDIR}/${chrom}/
	echo "#!/bin/bash" > $jobName
	echo "module load snpAD/0.3.10" >> $jobName
	echo "snpAD -c 12 -o ${OUTDIR}/${chrom}/${sample}.priors.txt -O ${OUTDIR}/${chrom}/${sample}.errors.txt $snpad" >> $jobName 
	chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 100G --time 24:00:00 -o ${out}/Rep30_${chrom}_${sample}.log --job-name ${sample}_$chrom -- $jobName

done <  <(head -n 1 repeatVCF3)
