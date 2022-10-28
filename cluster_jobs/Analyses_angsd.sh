#Obtain genotype likelihoods on all samples

module load angsd/0.939

# all sites with all samples
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/PCA/
jobName=$dir/out/pca_allsamples.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 28 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.0001 -SNP_pval 1e-6 -rf /home/qvw641/CottonTop_Tamarins/cluster_jobs/Chrom_autosomes_angsd  -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/CTT_allsamples" >> $jobName
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_all.log --job-name angsd -- $jobName

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

# PCA angsd
#Covariance matrix calculation
conda activate /projects/mjolnir1/apps/conda/pcangsd-0.98.2
pcangsd -beagle CTT_allsamples.beagle.gz -o CTT_allsamples_covmatrix -threads 10
pcangsd -beagle CTT_only.beagle.gz  -o CTT_only_covmatrix -threads 10
pcangsd -beagle CTT_only_maf.beagle.gz -o CTT_only_maf_covmatrix -threads 10

# Inbreeding with NgsRelate
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/
jobName=$dir/out/angsd_inbreeding.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_Bams_onlyCTT -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Inbreeding/CTT_only_rel" >> $jobName 
sbatch -c 1 --mem-per-cpu 50G --time 200:00:00 -o ${dir}/out/ANGSD_onlyCTT_rel.log --job-name angsdRel_CTT -- $jobName

zcat CTT_only_rel.mafs.gz | cut -f6 | sed 1d > CTT_rel.freq
/apps/NGSRELATE/2.0/GCC/bin/ngsRelate  -g CTT_rel.glf.gz -n 34 -f CTT_rel.freq > CTT_rel.ml


# Dissimilarity matrix ngsdist
# use beagle file produced for PCA
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/NgsDist/
mkdir -p $dir
mkdir -p ${dir}/out 
jobName=$dir/out/ngsdist_allsamples.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsDist/ngsDist --geno /projects/mjolnir1/people/qvw641/CottonTop/ANGSD/CTT_allsamples.beagle.gz  --n_ind 35 --n_sites 8702059 --probs TRUE --out ${dir}/Distance_allsamples.txt --n_boot_rep 3 --boot_block_size 1 --n_threads 5 --labels /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples" >> $jobName
sbatch -c 1 --mem-per-cpu 10G --time 200:00:00 -o ${dir}/out/NgsDist.log --job-name ngsdistCTT -- $jobName

for i in all;
do
	echo $i
	/projects/mjolnir1/apps/conda/fastme-2.1.6.1/bin/fastme -T 20 -i Distance_${i}.txt -s -D 3 -o ${i}.nwk
	head -n 1 ${i}.nwk > ${i}.main.nwk
	tail -n +2 ${i}.nwk | awk 'NF' > ${i}.boot.nwk
done

for i in all;
do      
        echo $i
        /projects/mjolnir1/apps/conda/raxml-ng-1.1.0/bin/raxml-ng  --support --tree ${i}.main.nwk --bs-trees ${i}.boot.nwk --prefix $i
done    

#Fst 
## per site
module load angsd/0.939
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/FST/Results_sites/
while read file;
do
        mkdir -p $dir
  	mkdir -p ${dir}/out/
	jobName=$dir/out/${file}.sh
        echo '#!/bin/bash' > $jobName
	echo "angsd -b /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_${file} -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta  \
                                  -out ${dir}/${file} -anc /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta \
                                 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                               -minMapQ 20 -minQ 20 -rf Chrom_autosomes -setMaxDepth 50 -doCounts 1 -GL 1 -doSaf 1" >> $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/Fst_${file}.log --job-name fst_${file} -- $jobName
done < ListSites


## per population TO DO
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/FST/Results_pop/
mkdir -p $dir
mkdir -p ${dir}/out/
while read file;
do
        jobName=$dir/out/${file}.sh
        echo '#!/bin/bash' > $jobName
        echo "angsd -b /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_${file} -ref /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta  \
                                  -out ${dir}/${file} -anc /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta \
                                 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                               -minMapQ 20 -minQ 20 -rf Chrom_autosomes -setMaxDepth 50 -doCounts 1 -GL 1 -doSaf 1" >> $jobName
        sbatch -c 1 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/Fst_${file}.log --job-name fst_${file} -- $jobName
done <  ListSitesPop




while read POP
do
        python3 ~/submit.py -c "/apps/ANGSD/0.931/bin/realSFS Results_sites/${POP}.saf.idx -P 8 > Results_sites/${POP}.sfs" -u 16 -n $POP -i . -w 8:00:00 -e out/$POP.err -o out/$POP.out

done < ListSites


# 2D-SFS between all sites 
module load angsd/0.939
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/FST/
mkdir -p ${dir}Results_sites_2D
while read  POP3
        do
                while read POP
                do
                echo $POP
                jobName=$dir/Results_sites_2D/out/2D_${POP3}_${POP}.sh
        	echo '#!/bin/bash' > $jobName
		echo "realSFS ${dir}Results_sites/${POP}.saf.idx ${dir}Results_sites/${POP3}.saf.idx -P 16  > ${dir}Results_sites_2D/$POP.$POP3.sfs" >> $jobName
	        sbatch -c 1 --mem-per-cpu 200G --time 8:00:00 -o ${dir}/Results_sites_2D/out/2D_${POP3}_${POP}.log --job-name 2s_${POP3}_${POP} -- $jobName
	done <  ListSites
done < ListSites

# Global Fst
dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/FST/
mkdir -p ${dir}Results_sites_FST
mkdir -p ${dir}Results_sites_FST/out/
while read  POP3
        do
        while read POP
               do
        	echo $POP
                jobName=$dir/Results_sites_FST/out/Fst_${POP3}_${POP}.sh
                echo '#!/bin/bash' > $jobName
	       	echo "realSFS fst index ${dir}Results_sites/${POP}.saf.idx ${dir}Results_sites/${POP3}.saf.idx -sfs ${dir}Results_sites_2D/${POP}.${POP3}.sfs -fstout ${dir}Results_sites_FST/$POP.$POP3" >>$jobName
	       sbatch -c 1 --mem-per-cpu 400G --time 12:00:00 -o ${dir}/Results_sites_FST/out/FST_${POP3}_${POP}.log --job-name FST_${POP3}_${POP} -- $jobName
        done < ListSites
done < ListSites

while read line;
do
        /apps/ANGSD/0.931/bin/realSFS fst stats Results_sites_FST/${line}.fst.idx > Results_sites_FST/${line}.fst
done < Results_sites_FST/ListComparison

while read file; do paste <(echo $file | cut -d"." -f1) <(echo $file | cut -d"." -f2) <(cat $file.fst); done < ListComparisons > Fst_total


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
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -setMaxDepth 50 -doCounts 1 \
                -GL 1 -doSaf 1" >> $jobName
	chmod +x $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/Het_${name}.log --job-name ha_${name} -- $jobName

done <  Samples_Bams


dir=/projects/mjolnir1/people/qvw641/CottonTop/ANGSD/Het/
while read sample;
do
	name=$(echo $sample | cut -d"/" -f 8 | cut -d"." -f1)
	jobName=$dir/out/SFS_${name}.sh
	echo '#!/bin/bash' > $jobName
	echo "realSFS ${dir}${name}.saf.idx -P 4 > ${dir}${name}.ml" >> $jobName
	chmod +x $jobName
	sbatch -c 1 --mem-per-cpu 600G --time 24:00:00 -o ${dir}/out/N2HetSFS_${name}.log --job-name sfs_${name} -- $jobName
done <  Samples_Bams
