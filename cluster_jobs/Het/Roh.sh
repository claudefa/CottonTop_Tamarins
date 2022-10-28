#!/bin/bash

indir='/projects/mjolnir1/people/qvw641/CottonTop/ROHs/'
outdir=${indir}'/'
errdir=${outdir}'/out'

mkdir -p $outdir $errdir 
file="/projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_onlyCTT_filter.vcf.gz" 

for i in 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e-9; 
do
	jobname=${errdir}/Bcftoolsroh_test_${i}.sh
	echo '#!/bin/bash' > $jobname
	#echo "bcftools roh -G30 --AF-dflt 0.0001  $file --rec-rate $i --skip-indels  -O r -o $outdir/BcfToolsRoh_alt0001recomb_${i}_final.txt" >> $jobname
	echo "bcftools roh -G30 --estimate-AF - $file --rec-rate $i --skip-indels  -O r -o $outdir/BcfToolsRoh_alt04recomb_${i}_final_AF.txt" >> $jobname
	echo $jobname
	chmod 755 $jobname
	sbatch -c 1 --mem-per-cpu 1G --time 3:00:00 -o ${errdir}/RohAF_${i}.log --job-name roh_${i} -- $jobname
done
