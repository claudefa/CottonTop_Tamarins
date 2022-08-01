#!/bin/bash
#Input Directories
dir='/scratch/devel/cfontser/CTT/VCF'
in_dir=$dir'/Filter'
#Output Directories
out_dir=$dir'/Filter'
err_dir=$out_dir'/out'
qu_dir=$out_dir'/qu'
mkdir -p $out_dir $err_dir $qu_dir
BCFTOOLS="/apps/BCFTOOLS/1.9/bin/bcftools"
vcfs=''
while read file;
do
        hello=$file
        vcfs=$vcfs" $hello"
done < <(ls ${in_dir}/CTT_*_filter.vcf.gz | sort -V)
echo $vcfs
jobname=${qu_dir}/concatsortVCFs.sh
echo "#!/bin/bash" > $jobname
echo "$BCFTOOLS concat $vcfs -Oz -o ${out_dir}/CTT_filter_old.vcf.gz; tabix -p vcf ${out_dir}/CTT_filter_old.vcf.gz;
	$BCFTOOLS reheader -h --samples Samples_reorder -Oz -o ${out_dir}/CTT_filter_final.vcf.gz ${out_dir}/CTT_filter_old.vcf.gz; tabix -p vcf ${out_dir}/CTT_filter_final.vcf.gz" >> $jobname
python3 ~/submit.py  -c $jobname -i . -e ${err_dir}/concatVCFs.err -o ${err_dir}/concatVCFs.out -n Concat -w 12:55:55 -u 1 
chmod 755 $jobname
