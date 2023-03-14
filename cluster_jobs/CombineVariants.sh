#!/bin/bash

#Ref genome
ref="/projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta"

#Directory
DIR="/projects/mjolnir1/people/qvw641/CottonTop/VCF/"
OUTDIR=$DIR"/MergeVCFs/"
out=${OUTDIR}"out/";
qu=${OUTDIR}"qu/";

# create directories
mkdir -p $OUTDIR
mkdir -p $out
mkdir -p $qu
module load bcftools

while read chrom;
do

	VCFs=""; while read line; do VCFs="$VCFs ${DIR}/$chrom/${line}_${chrom}.vcf.gz"; done  < <(cat Samples);
	jobName=${qu}/CTT_${chrom}.combine.sh
	echo "#!/bin/bash" > $jobName
	echo "bcftools merge --force-samples --merge all $VCFs -O z -o ${OUTDIR}/CTT_${chrom}.g.vcf.gz ;  tabix -p vcf ${OUTDIR}/CTT_${chrom}.g.vcf.gz" >> $jobName	
	chmod 755 $jobName
	echo $jobName
	sbatch -c 1 --mem-per-cpu 100G --time 24:00:00 -o ${out}/Comb_${chrom}.log --job-name Comb_$chrom -- $jobName

done < Chrom_autosomes
