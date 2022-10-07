#!/bin/bash
#Input Directories
dir='/projects/mjolnir1/people/qvw641/CottonTop/VCF'
in_dir=$dir'/FilteredVCF'
#Output Directories
out_dir=$dir'/Filter'
err_dir=$out_dir'/out'
qu_dir=$out_dir'/qu'
mkdir -p $out_dir $err_dir $qu_dir
module load bcftools
vcfs=''
while read file;
do
        hello=$file
        vcfs=$vcfs" $hello"
done < <(ls ${in_dir}/CTT_*_filter.vcf.gz | sort -V)
echo $vcfs
jobname=${qu_dir}/concatsortVCFs.sh
echo '#!/bin/bash' > $jobname
echo "module load bcftools; bcftools concat $vcfs -Oz -o ${out_dir}/CTT_allsamples_filter_nn.vcf.gz; tabix -p vcf ${out_dir}/CTT_allsamples_filter_nn.vcf.gz;
	$BCFTOOLS reheader -h --samples Samples -Oz -o ${out_dir}/CTT_allsamples_filter.vcf.gz ${out_dir}/CTT_allsamples_filter_nn.vcf.gz; tabix -p vcf ${out_dir}/CTT_allsamples_filter.vcf.gz" >> $jobname
chmod 755 $jobname

sbatch -c 1 --mem-per-cpu 100G --time 12:00:00 -o ${err_dir}/Concat.log --job-name Concat -- $jobname
