#Obtain genotype likelihoods on all samples

module load angsd/0.939

#Only CTT
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/PCA/
jobName=$dir/out/pca_CTTonlys.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref  /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 28 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.0001 -SNP_pval 1e-6 -rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd  -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/CTT_only" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT.log --job-name angsdCTT -- $jobName

# with Maf 0.05
#Only CTT
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/PCA/
jobName=$dir/out/pca_CTTonlysMaf.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref  /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 28 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -SNP_pval 1e-6 -rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd  -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/CTT_only_maf" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_maf.log --job-name angsdMCTT -- $jobName


# Only CTT historical
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/PCA/
jobName=$dir/out/pca_CTTonlysMaf_historical.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT_historical -ref  /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 21 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -SNP_pval 1e-6 -rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd  -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/CTT_only_maf_historical" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_maf.log --job-name angsdMCTT -- $jobName


# PCA angsd
#Covariance matrix calculation
conda activate /projects/mjolnir1/apps/conda/pcangsd-0.98.2
pcangsd -beagle CTT_only.beagle.gz  -o CTT_only_covmatrix -threads 10
pcangsd -beagle CTT_only_maf.beagle.gz -o CTT_only_maf_covmatrix -threads 10

# Inbreeding with NgsRelate
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel" >> $jobName 
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel.log --job-name angsdRel_CTT -- $jobName

zcat CTT_only_rel.mafs.gz | cut -f6 | sed 1d > CTT_rel.freq

dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/NGSrelate_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate -F 1  -g CTT_only_rel.glf.gz -n 34 -f CTT_rel.freq -O CTT_only_rel_inbreeding.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 48:00:00 -o ${dir}/out/NGSrelate_onlyCTT_rel.log --job-name Rel_CTT -- $jobName



# Inbreeding with only tv
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding_tvonly.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -noTrans 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel_noTrans" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel_notrans.log --job-name Rel_noTrans -- $jobName

zcat CTT_only_rel_noTrans.mafs.gz | cut -f6 | sed 1d > CTT_rel_notrans.freq

dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/NGSrelate_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate -F 1  -g CTT_only_rel.glf.gz -n 34 -f CTT_rel.freq -O CTT_only_rel_inbreeding.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 48:00:00 -o ${dir}/out/NGSrelate_onlyCTT_rel.log --job-name Rel_CTT -- $jobName

# Inbreeding with downsampled 5x
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding_down5x.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT_down5x -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel_5x" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel_5x.log --job-name Rel_5x -- $jobName

zcat CTT_only_rel_5x.mafs.gz | cut -f6 | sed 1d > CTT_rel_5x.freq
jobName=$dir/out/NGSrelate_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate -F 1  -g CTT_only_rel_5x.glf.gz -n 34 -f CTT_rel_5x.freq -O CTT_only_rel_inbreeding_5x.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 48:00:00 -o ${dir}/out/NGSrelate_onlyCTT_rel_5x.log --job-name Rel_CTT -- $jobName

# Inbreeding with downsampled 5x tv only
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding_tvonly_down5x.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT_down5x -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -noTrans 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel_5x_noTrans" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel_5x_notrans.log --job-name Rel_5x_noTrans -- $jobName


# Heterozygosity
# Calculate heterozygosity from SFS
module load angsd/0.939
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
while read sample;
do
	name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
	grep $name Samples_Bams > ${dir}${name}.tmp
	jobName=$dir/out/${name}.sh
	echo '#!/bin/bash' > $jobName
	echo  "angsd  -b ${dir}${name}.tmp -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta \
       		-rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd \
		-anc /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta  -out ${dir}${name} \
                -noTrans 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -setMaxDepth 50 -doCounts 1 \
                -GL 1 -doSaf 1" >> $jobName
	chmod +x $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/Het_${name}.log --job-name ha_${name} -- $jobName

done <  Samples_Bams_onlyCTT


dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
while read sample;
do
	name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
	jobName=$dir/out/SFS_${name}.sh
	echo '#!/bin/bash' > $jobName
	echo "realSFS ${dir}${name}.saf.idx -P 4 > ${dir}${name}.ml" >> $jobName
	chmod +x $jobName
	sbatch -c 1 --mem-per-cpu 600G --time 24:00:00 -o ${dir}/out/N2HetSFS_${name}.log --job-name sfs_${name} -- $jobName
done <  Samples_Bams_onlyCTT

# Fst
# see FST/fst.sh
