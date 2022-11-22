# Create vcf per scaffold and sample
./snpAD_1.sh
./snpAD_2.sh
# Combine variants per scaffold
./CombineVariants.sh 

# Filter VCF
./filtervcf.sh

# Concatenate all scaffolds
./concatVCFs_filter.sh

# Filter VCF to keep only CTT
jobName=/projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/out/filteronlyCTT.sh
echo '#!/bin/bash' > $jobName
echo "vcftools --gzvcf /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_allsamples_filter.vcf.gz --keep Samples_CTTonly --maf 0.000001 --minDP 3 --maxDP 50 --minGQ 30 --recode --recode-INFO-all --stdout | bgzip -c > /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_onlyCTT_filter.vcf.gz; tabix -p vcf /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_onlyCTT_filter.vcf.gz" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 6:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/out/Filter.log  --job-name filter -- $jobName


# VCF to plink - maf 0.05
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf CTT_allsamples_filter.vcf.gz --allow-extra-chr --double-id --maf 0.05 --geno 0.2 --recode --out CTT_allsamples_filter
sed -i 's/CM038391.1/1/g' CTT_allsamples_filter.map
sed -i 's/CM038392.1/2/g' CTT_allsamples_filter.map
sed -i 's/CM038393.1/3/g' CTT_allsamples_filter.map
sed -i 's/CM038394.1/4/g' CTT_allsamples_filter.map
sed -i 's/CM038395.1/5/g' CTT_allsamples_filter.map
sed -i 's/CM038396.1/6/g' CTT_allsamples_filter.map
sed -i 's/CM038397.1/7/g' CTT_allsamples_filter.map
sed -i 's/CM038399.1/8/g' CTT_allsamples_filter.map
sed -i 's/CM038400.1/9/g' CTT_allsamples_filter.map
sed -i 's/CM038401.1/10/g' CTT_allsamples_filter.map
sed -i 's/CM038402.1/11/g' CTT_allsamples_filter.map
sed -i 's/CM038403.1/12/g' CTT_allsamples_filter.map
sed -i 's/CM038404.1/13/g' CTT_allsamples_filter.map
sed -i 's/CM038405.1/14/g' CTT_allsamples_filter.map
sed -i 's/CM038406.1/15/g' CTT_allsamples_filter.map
sed -i 's/CM038407.1/16/g' CTT_allsamples_filter.map
sed -i 's/CM038408.1/17/g' CTT_allsamples_filter.map
sed -i 's/CM038409.1/18/g' CTT_allsamples_filter.map
sed -i 's/CM038410.1/19/g' CTT_allsamples_filter.map
sed -i 's/CM038411.1/20/g' CTT_allsamples_filter.map
sed -i 's/CM038412.1/21/g' CTT_allsamples_filter.map
sed -i 's/CM038413.1/22/g' CTT_allsamples_filter.map
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --allow-extra-chr --double-id --recode --make-bed --out CTT_allsamples_filter 

# VCF to plink - with singletons
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf CTT_allsamples_filter.vcf.gz --allow-extra-chr --double-id --maf 0.01 --geno 0.2 --recode --out CTT_allsamples_filter_singletons
sed -i 's/CM038391.1/1/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038392.1/2/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038393.1/3/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038394.1/4/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038395.1/5/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038396.1/6/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038397.1/7/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038399.1/8/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038400.1/9/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038401.1/10/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038402.1/11/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038403.1/12/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038404.1/13/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038405.1/14/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038406.1/15/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038407.1/16/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038408.1/17/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038409.1/18/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038410.1/19/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038411.1/20/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038412.1/21/g' CTT_allsamples_filter_singletons.map
sed -i 's/CM038413.1/22/g' CTT_allsamples_filter_singletons.map
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter_singletons --allow-extra-chr  --double-id --recode --make-bed --out CTT_allsamples_filter_singletons_LD

# Admixture only CTT
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_onlyCTT_filter.vcf.gz --maf 0.05 --geno 0.2 --allow-extra-chr --double-id --recode --out CTT_onlyCTT_filter
sed -i 's/CM038391.1/1/g' CTT_onlyCTT_filter.map
sed -i 's/CM038392.1/2/g' CTT_onlyCTT_filter.map
sed -i 's/CM038393.1/3/g' CTT_onlyCTT_filter.map
sed -i 's/CM038394.1/4/g' CTT_onlyCTT_filter.map
sed -i 's/CM038395.1/5/g' CTT_onlyCTT_filter.map
sed -i 's/CM038396.1/6/g' CTT_onlyCTT_filter.map
sed -i 's/CM038397.1/7/g' CTT_onlyCTT_filter.map
sed -i 's/CM038399.1/8/g' CTT_onlyCTT_filter.map
sed -i 's/CM038400.1/9/g' CTT_onlyCTT_filter.map
sed -i 's/CM038401.1/10/g' CTT_onlyCTT_filter.map
sed -i 's/CM038402.1/11/g' CTT_onlyCTT_filter.map
sed -i 's/CM038403.1/12/g' CTT_onlyCTT_filter.map
sed -i 's/CM038404.1/13/g' CTT_onlyCTT_filter.map
sed -i 's/CM038405.1/14/g' CTT_onlyCTT_filter.map
sed -i 's/CM038406.1/15/g' CTT_onlyCTT_filter.map
sed -i 's/CM038407.1/16/g' CTT_onlyCTT_filter.map
sed -i 's/CM038408.1/17/g' CTT_onlyCTT_filter.map
sed -i 's/CM038409.1/18/g' CTT_onlyCTT_filter.map
sed -i 's/CM038410.1/19/g' CTT_onlyCTT_filter.map
sed -i 's/CM038411.1/20/g' CTT_onlyCTT_filter.map
sed -i 's/CM038412.1/21/g' CTT_onlyCTT_filter.map
sed -i 's/CM038413.1/22/g' CTT_onlyCTT_filter.map
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_onlyCTT_filter --allow-extra-chr --double-id --recode --make-bed --out CTT_onlyCTT_filter

jobName=/projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/out/admixturesonlyCTT.sh
echo '#!/bin/bash' > $jobName
for j in {1..10};
do
for i in {1..8};
do
        echo "module load admixture; cd /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/; \
	admixture -s $j --cv /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/CTT_onlyCTT_filter.bed $i | tee /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/log_onlyCTT_${j}_${i}.out" >> $jobName
done
done
chmod +x $jobName
sbatch -c 1 --mem-per-cpu 2G --time 6:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/out/Adm_onlyCTT_%A_%a.log --array=1-80%10 --job-name admix -- $jobName


# Het and PCA
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --ibc --out CTT_allsamples_filter --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --missing --out CTT_allsamples_filter --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --pca --out CTT_allsamples_filter --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --het --out CTT_allsamples_filter --allow-extra-chr --double-id


/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter_singletons --ibc --out CTT_allsamples_filter_singletons --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter_singletons --test-missing --out CTT_allsamples_filter_singletons --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter_singletons --pca --out CTT_allsamples_filter_singletons --allow-extra-chr --double-id

# Basic stats and PCA only CTT
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly  --ibc --geno 0.2 --maf 0.05  --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly  --het --maf 0.05 --geno 0.2  --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly --maf 0.05 --missing --geno 0.2 --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly --maf 0.05 --pca --geno 0.2 --out CTT_filter_onlyCTT --allow-extra-chr --double-id

#EEMS
# apply eems to a set of historical cotton-top tamarins
EEMS/eems.sh

# F3 only CTT and with S.geoffroyi 
F3/f3_newsample.sh # with S. geoffroyi
F3/f3_onlyCTT.sh # only CTT


# Heterozygosity for all samples
# global scaffold
while read line; do sample=$(echo $line | awk '{print $1}'); id=$(echo $line | awk '{print $2}'); bash Het/calculateHet_scaffold_real.sh $sample $id; done < <(tail -n 15 Het/Samples)
while read line; do sample=$(echo $line | awk '{print $1}'); cat $sample/*het > ${sample}_total.het ; done < /home/qvw641/CottonTop_Tamarins/cluster_jobs/Het/Samples
# windows
Het/makewindows.sh
while read line; 
do  
	sample=$(echo $line | awk '{print $1}')
	id=$(echo $line | awk '{print $2}')
	jobName=/projects/mjolnir1/people/qvw641/CottonTop/VCF/Het/out/${sample}.window.sh
	echo '#!/bin/bash' >$jobName
	echo "bash Het/calculateHet_windows_real.sh 100000 $line" >> $jobName
	chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 1G --time 48:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/VCF/Het/out/WindHet_${sample}.log --job-name w_${sample}_$chrom -- $jobName
done < Het/Samples


# Allele imbalance
while read line; 
do
        sample=$(echo $line | awk '{print $1}')
	jobName=/projects/mjolnir1/people/qvw641/CottonTop/VCF/AB/out/${sample}.ab.sh
        echo '#!/bin/bash' > $jobName
	echo "vcftools --gzvcf /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/CTT_allsamples_filter.vcf.gz --indv $line --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | grep '0/1' | awk '{print \$10}' | cut -d':' -f2,3,4,5,6 | tr '[:]' '[\t]' > /projects/mjolnir1/people/qvw641/CottonTop/VCF/AB/AB_${line}.txt" >> $jobName
	chmod 755 $jobName
        sbatch -c 1 --mem-per-cpu 1G --time 2:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/VCF/AB/out/AB_${sample}.log --job-name ab_${sample} -- $jobName
done </home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples


# ROHs with BCFtools
# remove first the geoffroyi individual
./Het/Roh.sh

# ROHs with RoHan
./Rohan/rohan.sh

# Genetic Load
./counts_load.sh
