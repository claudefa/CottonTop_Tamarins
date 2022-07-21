# Obtain genotype likelihoods on the first scaffold

# all sites
python3 ~/submit.py -c "/apps/ANGSD/0.916/angsd -bam Samples -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 17 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.01 -SNP_pval 1e-6 -r CM038391.1: -out CTT_CM038391.1" -e angsd.err -o angsd.out -n angsd -t 10 -u 4 -r lowprio -w 50:00:00
python3 ~/submit.py -c "/apps/ANGSD/0.916/angsd -bam Samples -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 17 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.055555P_pval 1e-6 -r CM038391.1: -out CTT_CM038391.1_maf0.05" -e angsd.err -o angsd.out -n angsd -t 10 -u 4 -r lowprio -w 50:00:00

# without transversions
python3 ~/submit.py -c "/apps/ANGSD/0.916/angsd -bam Samples -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -rmTrans 1 -C 50 -baq 1 -minInd 17 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.01 -SNP_pval 1e-6 -r CM038391.1: -out CTT_CM038391.1_rmTrans" -e angsd_trans.err -o angsd_trans.out -n angsd_trans -t 10 -u 4 -r lowprio -w 50:00:00

# trimming first 10 bp
python3 ~/submit.py -c "/apps/ANGSD/0.916/angsd -bam Samples -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 10  -C 50 -baq 1 -minInd 17 -skipTriallelic 1 -GL 2 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMaf 0.01 -SNP_pval 1e-6 -r CM038391.1: -out CTT_CM038391.1_trim" -e angsd_trim.err -o angsd_trim.out -n angsd_trim -t 10 -u 4 -r lowprio -w 50:00:00

# PCA angsd
#Covariance matrix calculation
module unload PYTHON/3.4.3
module load PYTHON/2.7.14

/apps/PYTHON/2.7.14/AVX0/GCC/bin/python /apps/PCANGSD/20180209/pcangsd.py -beagle CTT_CM038391.1.beagle.gz -n 34 -o CTT_CM038391.1_covmatrix -threads 10
/apps/PYTHON/2.7.14/AVX0/GCC/bin/python /apps/PCANGSD/20180209/pcangsd.py -beagle CTT_CM038391.1_maf0.05.beagle.gz -n 34 -o CTT_CM038391.1_maf0.05_covmatrix -threads 10
/apps/PYTHON/2.7.14/AVX0/GCC/bin/python /apps/PCANGSD/20180209/pcangsd.py -beagle CTT_CM038391.1_rmTrans.beagle.gz -n 34 -o CTT_CM038391.1_rmTrans_covmatrix -threads 10
/apps/PYTHON/2.7.14/AVX0/GCC/bin/python /apps/PCANGSD/20180209/pcangsd.py -beagle CTT_CM038391.1_trim.beagle.gz -n 34 -o CTT_CM038391.1_trim_covmatrix -threads 10



#Relatedness
python3 ~/submit.py -c  "/apps/ANGSD/0.916/angsd -bam Samples -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta -uniqueOnly 1 -r CM038391.1: -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -out CTT_rel" -i . -e rel.err -o rel.out -n rel -u 2 -w 24:00:00

zcat CTT_rel.mafs.gz | cut -f6 | sed 1d > CTT_rel.freq
/apps/NGSRELATE/2.0/GCC/bin/ngsRelate  -g CTT_rel.glf.gz -n 34 -f CTT_rel.freq > CTT_rel.ml

#Fst
while read file;
do
        mkdir -p Results_sites
        python3 ~/submit.py -c "/apps/ANGSD/0.931/bin/angsd -b Samples_${file} -ref /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta  \
                                  -out Results_sites/${file} -anc /scratch/devel/cfontser/CTT/Assembly/SaguinusMidas_NCBI.fasta \
                                 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                               -minMapQ 20 -minQ 20 -rf Chrom_autosomes -setMaxDepth 50 -doCounts 1 \
                             -GL 1 -doSaf 1" -n $file -e out/${file}.err -o out/${file}.out -u 1 -r lowprio -w 48:00:00

done < ListSites

while read POP
do
        python3 ~/submit.py -c "/apps/ANGSD/0.931/bin/realSFS Results_sites/${POP}.saf.idx -P 8 > Results_sites/${POP}.sfs" -u 16 -n $POP -i . -w 8:00:00 -e out/$POP.err -o out/$POP.out

done < ListSites


# 2D-SFS between all sites only 1st scaffold

while read  POP3
        do
                while read POP
                do
                echo $POP
                mkdir -p Results_sites_2D
                python3 ~/submit.py -c "/apps/ANGSD/0.931/bin/realSFS Results_sites/${POP}.saf.idx Results_sites/${POP3}.saf.idx -P 16  > Results_sites_2D/$POP.$POP3.sfs" -u 20 -t 2 -n $POP -i . -w 12:00:00 -e out/$POP.$POP3.err -o out/$POP.$POP3.out
        done < ListSites
done < <(head -n 1 ListSites)

# Global Fst

while read  POP3
        do
        while read POP
               do
                mkdir -p Results_sites_FST
                python3 ~/submit.py -c "/apps/ANGSD/0.931/bin/realSFS fst index Results_sites/${POP}.saf.idx Results_sites/${POP3}.saf.idx -r CM038391.1 -sfs Results_sites_2D/${POP}.${POP3}.sfs -fstout Results_sites_FST/$POP.$POP3" -u 1  -n ${POP}_${POP3} -i . -w 8:00:00 -e out/$POP.$POP3.FST.err -o out/$POP.$POP3.FST.out
        done < ListSites
done < ListSites

while read line;
do
        /apps/ANGSD/0.931/bin/realSFS fst stats Results_sites_FST/${line}.fst.idx > Results_sites_FST/${line}.fst
done < Results_sites_FST/ListComparison

while read file; do paste <(echo $file | cut -d"." -f1) <(echo $file | cut -d"." -f2) <(cat $file.fst); done < ListComparisons > Fst_total



