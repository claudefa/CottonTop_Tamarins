#!/bin/bash

in="/projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/"
vcf=${in}"CTT_allsamples_filter.vcf.gz"
outdir="/projects/mjolnir1/people/qvw641/CottonTop/F3/"
qu="/projects/mjolnir1/people/qvw641/CottonTop/F3/qu/"
out="/projects/mjolnir1/people/qvw641/CottonTop/F3/out/"
script="/home/qvw641/CottonTop_Tamarins/cluster_jobs/F3/"

mkdir -p $outdir
mkdir -p $qu
mkdir -p $out

module load bcftools
jobName1=$qu/eigensoftNewVCF.lst
echo "#!/bin/bash" > $jobName1
echo "bcftools view -i 'F_MISSING<0.2 && MAF > 0.00001'  $vcf | bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' | python ${script}vcf_to_eigensoft.py > ${outdir}CTT_newvcf.eigensoft" >> $jobName1
chmod 755 $jobName1
#sbatch -c 1 --mem-per-cpu 100G --time 10:00:00 -o ${out}/EigenNV.log --job-name eigen -- $jobName1

#manually add a 0 column at the end to account for the reference: G. midas

jobName2=$qu/snpNewVCF.lst
echo '#!/bin/bash' >$jobName2
echo "bcftools view -i 'F_MISSING<0.2 && MAF > 0.00001'  $vcf |  awk '/^#/{next}{print \"\t\"\$1\":\"\$2\"\t1\t0\t\"\$2\"\t\"\$4\"\t\"\$5 }' > ${outdir}CTT_newvcf.snp" >> $jobName2
chmod 755 $jobName2
#sbatch -c 1 --mem-per-cpu 100G --time 4:00:00 -o ${out}/SnpNV.log --job-name snp -- $jobName2

# obtain the comparisons
#awk '{print $3}' ${outdir}CTT_newvcf.ind | sort | uniq > ${outdir}sitesNV
#while read i ; do while read j; do echo  $i $j "Reference";done < ${outdir}sitesNV; done < ${outdir}sitesNV > ${outdir}list_sitesNV_qp3

module load admixtools/
jobName3=$qu/f3_commandNewVCF.sh
echo '#!/bin/bash' >$jobName3
echo 'module load admixtools/' >> $jobName3
echo "qp3Pop -p ${scripts}CTT_newVCF.par > ${outdir}One_logfile_Sites_CTT_newvcf_repeat2" >> $jobName3  
chmod 755 $jobName3
sbatch -c 1 --mem-per-cpu 10G --time 1:00:00 -o ${out}/Parf3.log --job-name f3_newVCF -- $jobName3

#grep 'result:' ${outdir}One_logfile_Sites_CTT_newvcf* | awk '{print $3, $4, $5, $6, $7, $8, $9, $10}' | sort | uniq > ${outdir}all_sites_newvcf.qp3Pop.out
