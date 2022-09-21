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

while read chrom;
do
# send array
array=$qu/SnpAD_${chrom}_array.sh

echo "#!/bin/bash" > $array
echo "#" >> $array
echo "#SBATCH -c 8" >> $array
echo "#SBATCH --mem-per-cpu 200G" >> $array
echo "#SBATCH --time=30:00:00"  >> $array
echo "#SBATCH --array=1-35%10" >> $array
echo "#SBATCH --output=${out}/snpad.${chrom}.%A_%a.log" >> $array
echo "#SBATCH --job-name snpad.${chrom}" >> $array

echo "conda activate /projects/mjolnir1/apps/conda/picard-2.15.0/" >> $array
echo "sample=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p ${DIR}/Samples | awk '{print \$1}')" >> $array
echo "echo \$sample" >> $array
echo "bam=${IN}\${sample}.SagMidas.realigned.bam" >> $array
echo "snpad=${OUTDIR}/${chrom}/\${sample}_${chrom}.snpAD" >> $array
echo "mkdir -p ${OUTDIR}/${chrom}/" >> $array
echo "module load snpAD/0.3.10" >> $array
echo "Bam2snpAD -r $chrom -f $ASSEMBLY -Q 30 \$bam  > \$snpad; \
      snpAD -c 12 -o ${OUTDIR}/${chrom}/\${sample}.priors.txt -O ${OUTDIR}/${chrom}/\${sample}.errors.txt \$snpad ; \
      snpADCall -N 1 -e ${OUTDIR}/${chrom}/\${sample}.errors.txt -p \"\`cat ${OUTDIR}/${chrom}/\${sample}.priors.txt\`\" \$snpad > \$vcf; \
      bgzip \$vcf; tabix -p vcf \${vcf}.gz; rm \$snpad; rm ${OUTDIR}/${chrom}/\${sample}.priors.txt; rm ${OUTDIR}/${chrom}/\${sample}.errors.txt" >> ${array}


#sbatch $array
done < Chrom_autosomes
