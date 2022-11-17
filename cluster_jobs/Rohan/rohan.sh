source /home/qvw641/bin/rohan_env/bin/activate

# run rohan
dir=/projects/mjolnir1/people/qvw641/CottonTop/ROHan

while read sample;
do 
	jobName=${dir}/out/jobName_${sample}.sh
	echo '#!/bin/bash' > $jobName
	echo "/home/qvw641/bin/rohan/src/rohan -t 16 --tstv 2.02 --size 1000000 --rohmu 2.5e-4 --auto Chrom_autosomes -o /projects/mjolnir1/people/qvw641/CottonTop/ROHan/${sample}_2.5e4_allscaffolds_1Mb /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta /projects/mjolnir1/people/qvw641/CottonTop/BAMs/${sample}.SagMidas.autosome.bam" >> $jobName
	cat $jobName | sbatch -c 1 --mem-per-cpu 1G --time 72:00:00  -o ${dir}/out/2e4_1Mb_Rohan_$sample.log --job-name e4_$sample.Rohan --            

done <  Samples
