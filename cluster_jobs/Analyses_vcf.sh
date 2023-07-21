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
for i in {9..10};
do
        echo "module load admixture; cd /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/; \
	admixture -s $j --cv /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/CTT_onlyCTT_filter.bed $i | tee /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/Rep${j}/log_onlyCTT_${j}_${i}.out" >> $jobName
done
done
chmod +x $jobName
sbatch -c 1 --mem-per-cpu 2G --time 6:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter/admixture/out/Adm_onlyCTT_%A_%a.log --array=1-20%10 --job-name admix -- $jobName


# Basic stats and PCA only CTT
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly  --ibc --geno 0.2 --maf 0.05  --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly  --het --maf 0.05 --geno 0.2  --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly --maf 0.05 --missing --geno 0.2 --out CTT_filter_onlyCTT --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --keep /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_CTTonly --maf 0.05 --pca --geno 0.2 --out CTT_filter_onlyCTT --allow-extra-chr --double-id

#EEMS
# apply eems to a set of historical cotton-top tamarins
EEMS/eems.sh

# F3  
F3/f3_onlyCTT.sh 

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

# ROHs with RoHan
./rohan.sh

# Genetic Load
java -Xmx8g -jar /projects/mjolnir1/apps/snpeff-5.1d/snpEff.jar -c snpEff.config -v SaguinusMidas_NCBI rename.vcf.gz > CTT_allsamples.ann.vcf
vcftools --vcf CTT_allsamples.ann.vcf --keep samples.list --maf 0.0001 --recode --recode-INFO-all --stdout | bgzip -c > CTT_oedipusSamples.ann.vcf.gz
vcftools --gzvcf CTT_oedipusSamples.ann.vcf.gz --keep Samples_5x --max-missing 1 --maf 0.0000001 --stdout --recode --recode-INFO-all | bgzip -c > CTT_5x.ann.vcf.gz;
./GeneticLoad/counts_load_dp5.sh

# Distribution of read depth and genotype quality per sample 
while read line;
do
        sample=$(echo $line | awk '{print $1}')
        jobName=/projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/out/${sample}.quality.sh
        echo '#!/bin/bash' > $jobName
        echo "vcftools --gzvcf /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/CTT_5x.ann.vcf.gz --indv $sample --minDP 5 --maxDP 50 --minGQ 30 --stdout --recode --recode-INFO-all | grep -v \"#\" |  grep HIGH | awk '{print \$10}' | cut -d':' -f1,2,8 | tr '[,]' '[\t]'  > /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/DP_GQ_${sample}_HIGH_v2.txt" >>  $jobName
        chmod 755 $jobName
        sbatch -c 1 --mem-per-cpu 1G --time 8:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/out/DP_GQ_${sample}_high.log --job-name gq_${sample} -- $jobName
        echo "vcftools --gzvcf /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/CTT_5x.ann.vcf.gz --indv $sample --minDP 5 --maxDP 50 --minGQ 30 --stdout --recode --recode-INFO-all | grep -v \"#\" |  grep MODERATE | awk '{print \$10}' | cut -d':' -f1,2,8 | tr '[,]' '[\t]'  > /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/DP_GQ_${sample}_MODERATE_v2.txt" >>  $jobName
        chmod 755 $jobName
        sbatch -c 1 --mem-per-cpu 1G --time 8:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/out/DP_GQ_${sample}_moderate.log --job-name gq_${sample} -- $jobName
        echo "vcftools --gzvcf /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/CTT_5x.ann.vcf.gz --indv $sample --minDP 5 --maxDP 50 --minGQ 30 --stdout --recode --recode-INFO-all | grep -v \"#\" |  grep synonymous_variant | awk '{print \$10}' | cut -d':' -f1,2,8 | tr '[,]' '[\t]'  > /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/DP_GQ_${sample}_SYNONYM_v2.txt" >>  $jobName
        chmod 755 $jobName
        sbatch -c 1 --mem-per-cpu 1G --time 8:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/GeneticLoad/out/DP_GQ_${sample}_synonym_v2.log --job-name gq_${sample} -- $jobName
done < /home/qvw641/CottonTop_Tamarins/cluster_jobs/Samples_5x

