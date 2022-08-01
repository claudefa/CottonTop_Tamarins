#!/bin/bash

sample=$1;
if [ -z "$sample" ];then printf "### Calculate Het by scaffold ###\n\nUSAGE:\n\nbash calculateHet_scaffold.sh Sample  \n\n";exit;fi

dir="/scratch/devel/cfontser/CTT/VCF/Heterozygosity"
het=${dir}/Het_real/${sample}/
out=${dir}/Het_real/out/
qu=${dir}/Het_real/qu/

mkdir -p $het $qu $out

while read line;
do 

	region=$(echo $line | awk '{print $1":"$2"-"$3}');
	read chrom start end <<<$(echo $line);
        jobName=${qu}/Scaffold_${chrom}_${sample}.sh
	vcf="/scratch/devel/cfontser/CTT/VCF/MergeVCFs/CTT_${chrom}.g.vcf.gz"

	echo "#!/bin/bash" > ${jobName} 
	echo "callable=\$(tabix -h $vcf $region | vcftools --vcf - --indv $sample --recode --max-missing 1 --recode-INFO-all --stdout --minDP 3 --maxDP 30 --minQ 30 | grep -v '#' | wc -l ) " >> ${jobName}
	echo "if [ \$callable -gt 0 ];">> ${jobName}
	echo "then" >> ${jobName}
	echo "	tabix -h $vcf $region | vcftools --vcf - --indv $sample --recode --recode-INFO-all --max-missing 1 --stdout --minDP 3 --maxDP 30 --minQ 30 | vcfhetcount | tail -n 1 | awk -v l=\$callable -v chrom=$chrom '{print chrom\"\t\"\$1\"\t\"l\"\t\"\$1/l}' > ${het}/${sample}_${chrom}.het" >> ${jobName}
	echo "else" >> ${jobName}
	echo "	echo $line | awk -v l=0 -v het=NA '{print \$1\"\t\"l\"\t\"l\"\t\"het}' > ${het}/${sample}_${chrom}.het" >> ${jobName}
	echo "fi" >> ${jobName}
	echo "callable=''" >> ${jobName} 
	chmod 755 $jobName
	python3 ~/submit.py  -u 1 -c $jobName -n Het_${sample} -o ${out}/Scaff_${sample}.out -e ${out}/Scaff${sample}.err -w "10:00:00"

done < <(tail -n 21 Windows/Assembly_scaffolds_autosomes.bed)
