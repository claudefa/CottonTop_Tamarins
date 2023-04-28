#Obtain genotype likelihoods on all samples

module load angsd/0.939

# Inbreeding with NgsRelate
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding_tvonly.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -noTrans 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.00001 -SNP_pval 1e-6 -P 4 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel_noTrans" >> $jobName
sbatch -c 4 --mem-per-cpu 10G --time 100:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel_notrans.log --job-name Rel_noTrans -- $jobName

zcat CTT_only_rel_noTrans.mafs.gz | cut -f6 | sed 1d > CTT_rel_notrans.freq

dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/NGSrelate_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate -F 1  -g CTT_only_rel_noTrans.glf.gz -n 34 -f CTT_rel_notrans.freq -O CTT_only_rel_notrans_inbreeding.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 12:00:00 -o ${dir}/out/NGSrelate_onlyCTT_rel.log --job-name Rel_CTT -- $jobName

# Inbreeding with downsampled 5x
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding_tvonly_down5x.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT_down5x -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -noTrans 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.00001 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel_5x_noTrans" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 48:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel_5x_notrans.log --job-name Rel_5x_noTrans -- $jobName

zcat CTT_only_rel_5x_noTrans.mafs.gz | cut -f6 | sed 1d > CTT_rel_5x_noTrans.freq
jobName=$dir/out/NGSrelate_inbreeding_tvonly_down5x.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate -F 1  -g CTT_only_rel_5x_noTrans.glf.gz -n 34 -f CTT_rel_5x_noTrans.freq -O CTT_only_rel_inbreeding_5x_noTrans.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 7:00:00 -o ${dir}/out/NGSrelate_onlyCTT_rel_5x_down.log --job-name Rel_CTTDown -- $jobName

# Heterozygosity
# Calculate heterozygosity from SFS
module load angsd/0.939
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
echo '#!/bin/bash' > $jobName
while read sample;
do
	name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
	grep $name Samples_Bams > ${dir}${name}.tmp
	jobName=$dir/out/${name}.sh
	echo '#!/bin/bash' > $jobName
	echo  "angsd  -b ${dir}${name}.tmp -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta \
       		-rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd \
		-anc /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta  -out ${dir}${name} \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -setMaxDepth 50 -doCounts 1 \
                -GL 1 -doSaf 1" >> $jobName
	chmod +x $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/Het_${name}.log --job-name ha_${name} -- $jobName

done <  Samples_Bams_onlyCTT

#with winSFS
module load winsfs/0.7.0
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
jobName=$dir/out/winsfs.sh
echo '#!/bin/bash' > $jobName
while read sample;
do
        name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
        echo "winsfs ${dir}${name}.saf.idx > ${dir}${name}_win3.ml" >> $jobName
done <  Samples_Bams_onlyCTT

cat $jobName | sbatch -c 1 --array=1 --mem-per-cpu 600G --time 8:00:00 -o ${dir}/out/winSFS_%A_%a.log --job-name winsfsT -- $jobName

# Calculate heterozygosity from SFS but with high coverage bam files downsampled to 5x
module load angsd/0.939
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
jobName=$dir/out/down_tv.sh
echo '#!/bin/bash' > $jobName
while read sample;
do
        name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
        grep $name Samples_Bams_onlyCTT_down5x > ${dir}${name}.tmp
        echo  "angsd  -b ${dir}${name}.tmp -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta \
                -rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd \
                -anc /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta  -out ${dir}Down/${name}_down_TV \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -setMaxDepth 50 -doCounts 1 \
                -GL 1 -doSaf 1" >> $jobName
done <  Samples_Bams_onlyCTT_down5x

cat $jobName | sbatch -c 1 --array=1-34%10 --mem-per-cpu 50G --time 80:00:00 -o ${dir}/out/Het_downTV_%A_%a.log --job-name ha_down_TV --

#with WinSFS
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
jobName=$dir/out/down_tv_winsfs.sh
echo '#!/bin/bash' > $jobName
while read sample;
do
        name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
        echo "winsfs ${dir}Down/${name}_down_TV.saf.idx > ${dir}Down/${name}_win.ml" >> $jobName
done <  Samples_Bams_onlyCTT_down5x

cat $jobName | sbatch -c 1 --array=1-34%10 --mem-per-cpu 1G --time 24:00:00 -o ${dir}/out/N2HetwinSFS_%A_%a.log --job-name winsfs -- $jobName

# Fst
# see FST/fst.sh
