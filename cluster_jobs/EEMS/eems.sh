/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf /projects/mjolnir1/people/vbz170/projects/CTT/Shotgun_paper/CTT_filter.vcf.gz --keep keep_historicals.txt --recode --out /projects/mjolnir1/people/qvw641/CottonTop/HistoricalCTT/EEMS/ --allow-extra-chr --double-id --maf 0.05 --geno 0.2
sed -i 's/CM038391.1/chr1/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038392.1/chr2/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038393.1/chr3/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038394.1/chr4/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038395.1/chr5/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038396.1/chr6/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038397.1/chr7/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038399.1/chr8/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038400.1/chr9/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038401.1/chr10/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038402.1/chr11/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038403.1/chr12/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038404.1/chr13/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038405.1/chr14/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038406.1/chr15/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038407.1/chr16/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038408.1/chr17/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038409.1/chr18/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038410.1/chr19/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038411.1/chr20/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038412.1/chr21/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map
sed -i 's/CM038413.1/chr22/g' /projects/mjolnir1/people/qvw641/CottonTop/EEMS/HistoricalCTT.map

/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file HistoricalCTT --make-bed --out HistoricalCTT
bed2diffs_v1 --bfile HistoricalCTT
module load gcc
echo '#!/bin/bash' > job1.txt
echo 'module load gcc' >> job1.txt
echo '/projects/mjolnir1/apps/bin/runeems_snps --params /home/qvw641/CottonTop_Tamarins/cluster_jobs/EEMS/params-chain1.ini' >> job1.txt
cat job1.txt | sbatch -c 1 --mem-per-cpu 100G --time 200:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/EEMS/out/EEMS1.log --job-name EEMS1

for i in {2..10};
do
	echo "#!/bin/bash" > job2.txt
	echo "module load gcc" >> job2.txt
	echo "/projects/mjolnir1/apps/bin/runeems_snps --params /home/qvw641/CottonTop_Tamarins/cluster_jobs/EEMS/params-chain${i}.ini" >> job2.txt
	cat job2.txt | sbatch -c 1 --mem-per-cpu 100G --time 200:00:00 -o /projects/mjolnir1/people/qvw641/CottonTop/EEMS/out/EEMS${i}.log --job-name EEMS${i}
done


