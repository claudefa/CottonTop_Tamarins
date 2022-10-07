# Create vcf per scaffold and sample
./snpAD.sh
./snpAD_2.sh
./snpAD_3.sh
# Combine variants per scaffold
./CombineVariants.sh 

# Filter VCF
./filtervcf.sh

# Concatenate all scaffolds
./concatVCFs_filter.sh

# Basic stats and PCA

vcftools --gzvcf CTT_filter.vcf.gz --het --maf 0.05 --out CTT_filter
plink --vcf CTT_filter.vcf.gz --ibc --maf 0.05  --out CTT_filter --allow-extra-chr --double-id
plink --vcf CTT_filter.vcf.gz --maf 0.05 --test-missing --out CTT_filter --allow-extra-chr --double-id
plink --vcf CTT_filter.vcf.gz --maf 0.05 --pca --out CTT_filter --allow-extra-chr --double-id

# projected PCA 
 plink --vcf CTT_CM038391.1_filter.vcf.gz --maf 0.05 --pca --out CTT_CM038391.1 --allow-extra-chr --pca-cluster-names Historical --within clusters

# admixture
plink --vcf CTT_CM038391.1_filter.vcf.gz  --recode12 --indep-pairwise 50 10 0.1 --maf 0.05 --out CTT_CM038391.1
for K in 1 2 3 4 5; do /apps/ADMIXTURE/1.23/bin/admixture --cv CTT_CM038391.1.ped $K | tee log${K}.out; done


#EEMS
# apply eems to a set of historical cotton-top tamarins
module load gcc/6.3.0
module load boost/latest
module load eems

python3 ~/submit.py -c "vcftools --gzvcf  /scratch/devel/cfontser/CTT/VCF/Filter/CTT_CM038391.1_filter.vcf.gz --keep historical.txt --maf 0.05 --minDP 3 --minGQ 20 --max-missing 0.8 --stdout --recode --recode-INFO-all | bgzip -c > HistoricalCTT_CM038391.1_filter.vcf.gz ;tabix -p vcf HistoricalCTT_CM038391.1_filter.vcf.gz" -i . -e filter.err -o filter.out -n filter -u 1 -w 2:00:00

plink --vcf HistoricalCTT_CM038391.1_filter.vcf.gz --recode --out HistoricalCTT --allow-extra-chr
sed -i 's/CM038391.1/chr1/g' HistoricalCTT.map
plink --file HistoricalCTT --make-bed --out HistoricalCTT

bed2diffs_v1 --bfile HistoricalCTT

python3 ~/submit.py -i . -e eems.err -o eems.out -u 1 -w 23:00:00 -n eems -c "runeems_snps --params params-chain1.ini"


# Heterozygosity
# global scaffold
while read line; do bash calculateHet_scaffold_real.sh $line; done < Samples

# windows
./makewindows.sh
while read line; do python3 ~/submit.py -c "bash calculateHet_windows_real.sh 100000 $line" -i . -e out/${line}_window.err -o out/${line}_window.out -n w$line -u 1 -w 6:00:00; done < Samples



# New VCF

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

# Admixture in /projects/mjolnir1/people/qvw641/CottonTop/VCF/Filter
for i in {2..10};
do
	admixture --cv CTT_allsamples_filter.bed $i | tee log${i}.out
done 

grep -h CV log*.out

# Het and PCA
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --ibc --out CTT_allsamples_filter --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --test-missing --out CTT_allsamples_filter --allow-extra-chr --double-id
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file CTT_allsamples_filter --pca --out CTT_allsamples_filter --allow-extra-chr --double-id



	
