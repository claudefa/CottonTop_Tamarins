# Create vcf per scaffold and sample
./snpAD.sh

# Combine variants per scaffold
./CombineVariants.sh 

# Filter VCF

while read line
do
python3 ~/submit.py -c "bcftools filter -e'FORMAT/DP<3 && FORMAT/GQ<30 && FORMAT/DP>50' /scratch/devel/cfontser/CTT/VCF/MergeVCFs/CTT_${line}.g.vcf.gz | bcftools view -m2 -M2 -v snps -i 'F_MISSING<0.3' | bgzip -c > CTT_${line}_filter.vcf.gz; tabix -p vcf CTT_${line}_filter.vcf.gz"  -i. -e out/filter.err -o out/filter.out -u 2 -w 12:00:00  -n filter
done < Chrom_autosomes

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

