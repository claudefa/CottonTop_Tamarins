#!/bin/bash

sample=$1;
id=$2
if [ -z "$sample" ];then printf "### Calculate Het by scaffold ###\n\nUSAGE:\n\nbash calculateHet_scaffold.sh Sample ID_vcf \n\n";exit;fi

dir=/projects/mjolnir1/people/qvw641/CottonTop/VCF/Het/
het=${dir}/Het_real/${sample}/
out=${dir}/Het_real/out/
qu=${dir}/Het_real/qu/

mkdir -p $het $qu $out

while read line;
do 

	region=$(echo $line | awk '{print $1":"$2"-"$3}');
	read chrom start end <<<$(echo $line);
        jobName=${qu}/Scaffold_${chrom}_${sample}.sh
	vcf="/projects/mjolnir1/people/qvw641/CottonTop/VCF/MergeVCFs/CTT_${chrom}.g.vcf.gz"

	echo "#!/bin/bash" > ${jobName} 
	echo "callable=\$(tabix -h $vcf $region | vcftools --vcf - --indv $id --recode --max-missing 1 --recode-INFO-all --stdout --minDP 3 --maxDP 30 --minQ 30 | grep -v '#' | wc -l ) " >> ${jobName}
	echo "if [ \$callable -gt 0 ];">> ${jobName}
	echo "then" >> ${jobName}
	echo "	tabix -h $vcf $region | vcftools --vcf - --indv $id --recode --recode-INFO-all --max-missing 1 --stdout --minDP 3 --maxDP 30 --minQ 30 | grep 0/1 | wc -l | awk -v l=\$callable -v chrom=$chrom '{print chrom\"\t\"\$1\"\t\"l\"\t\"\$1/l}' > ${het}/${sample}_${chrom}.het" >> ${jobName}
	echo "else" >> ${jobName}
	echo "	echo $line | awk -v l=0 -v het=NA '{print \$1\"\t\"l\"\t\"l\"\t\"het}' > ${het}/${sample}_${chrom}.het" >> ${jobName}
	echo "fi" >> ${jobName}
	echo "callable=''" >> ${jobName} 
	chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 1G --time 8:00:00 -o ${out}/Het_${chrom}_${sample}.log --job-name Het_${sample}_$chrom -- $jobName

done < <(tail -n 21 /home/qvw641/CottonTop_Tamarins/cluster_jobs/Het/Windows/Assembly_scaffolds_autosomes.bed)
