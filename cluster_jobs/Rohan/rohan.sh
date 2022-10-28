source /home/qvw641/bin/rohan_env/bin/activate

# run rohan
dir=/projects/mjolnir1/people/qvw641/CottonTop/ROHan

while read sample;
do 
	jobName=${dir}/out/jobName_${sample}.sh
	echo '#!/bin/bash' > $jobName
	echo "/home/qvw641/bin/rohan/src/rohan -t 8 --rohmu 2e-5 -o /projects/mjolnir1/people/qvw641/CottonTop/ROHan/${sample} /projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta /projects/mjolnir1/people/qvw641/CottonTop/BAMs/${sample}.SagMidas.autosome.bam --size 100000" >> $jobName
	cat $jobName | sbatch -c 1 --mem-per-cpu 10G --time 100:00:00  -o ${dir}/out/Rohan_$sample.log --job-name Rohan --            

done < <(head -n 10 ../Samples)
