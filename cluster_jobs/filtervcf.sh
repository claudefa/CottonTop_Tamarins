#!/bin/bash

#Script to filter variants per chromosomes 

#Directory
DIR="/projects/mjolnir1/people/qvw641/CottonTop/VCF/"
INDIR=$DIR"/MergeVCFs/"
OUTDIR=$DIR"/FilteredVCF/"
out=${OUTDIR}"out/";
qu=${OUTDIR}"qu/";

# create directories
mkdir -p $OUTDIR
mkdir -p $out
mkdir -p $qu


module load bcftools
while read line
do
	echo $line
	jobName=$qu/Filter_${line}.sh
        
	echo '#!/bin/bash' > $jobName
        echo "module load bcftools; bcftools filter -e 'FORMAT/DP<3 && FORMAT/GQ<30 && FORMAT/DP>50' ${INDIR}/CTT_${line}.g.vcf.gz | bcftools view -m2 -M2 -v snps -i 'F_MISSING<0.3' | bgzip -c > ${OUTDIR}/CTT_${line}_filter.vcf.gz; tabix -p vcf ${OUTDIR}/CTT_${line}_filter.vcf.gz" >> $jobName 
	chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 100G --time 12:00:00 -o ${out}/Filt_${line}.log --job-name filt$line -- $jobName	
done < <(tail -n 21 Chrom_autosomes)
