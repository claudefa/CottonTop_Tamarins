#!/bin/bash

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

#REF
ASSEMBLY="/projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta"
module load snpAD/0.3.10

while read sample
do
	while read chrom;
	do
        jobName=$qu/SnpAD_${chrom}_$sample.sh
        bam=${IN}${sample}.SagMidas.autosome.bam
        snpad=${OUTDIR}/${chrom}/${sample}_${chrom}.snpAD
        vcf=${OUTDIR}/${chrom}/${sample}_${chrom}.vcf
        mkdir -p ${OUTDIR}/${chrom}/
        echo "#!/bin/bash" > $jobName
        echo "module load snpAD/0.3.10" >> $jobName
        echo "snpADCall -N 1 -e ${OUTDIR}/${chrom}/${sample}.errors.txt -p \"\`cat ${OUTDIR}/${chrom}/${sample}.priors.txt\`\" $snpad > $vcf; \
      		bgzip $vcf; tabix -p vcf ${vcf}.gz;" >> $jobName
        chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 100G --time 2:00:00 -o ${out}/snpAD_f_${sample}_${chrom}.log --job-name ${chrom}_$sample -- $jobName

done < Chrom_autosomes	
done < Samples
