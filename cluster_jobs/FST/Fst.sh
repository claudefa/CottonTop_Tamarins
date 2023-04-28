#Fst with ANGSD/v.921

################################################ By Region Fst ######################################################################################

#Number of bams in each population:
# 9 north_HisPop.txt (San Juan, Coloso, Cauca) 
# 6 north_ModPop.txt (San Juan, Ciebal)
# 11 central_HisPop.txt (Tierralta, Planeta Rica, Caracas, Unknown)
# 6 south_HisPop.txt (Mutata, Turbo, Ajorna)
# 2 south_ModPop.txt (Tulenpaa)

#make a bash file per population with the following
#run with  sbatch -c 5 --mem-per-cpu 4G --
# -P is number of threads
# 1. Estimate GLs

cd $STATS/Fst/Fst_region

/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam south_HisPop.txt -doSaf 1 -out south_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 3

/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam south_ModPop.txt -doSaf 1 -out south_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 2


# 2. Estimate 2-dimensional frequency spectrum from the site allele frequency likelihoods (for each pair, 3 in this case)
#create bash files with:

# Historical South vs Modern South
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS south_HisPop.saf.idx south_ModPop.saf.idx >HisSouth_ModSouth.ml

# run the above with sbatch
sbatch -c 4 --mem-per-cpu 20G --time=7-00 -J HisCentral_ModNorth_SFS  -- scripts/runSFS_HisCentral_ModNorth.sh


# 3. pops fsts - send all below through sbtach --> sbatch -c 1 --mem-per-cpu 4G --time=7-00 -J



/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index south_HisPop.saf.idx south_ModPop.saf.idx  -fstout South_only -whichFst 1  -sfs HisSouth_ModSouth.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats South_only.fst.idx


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index south_HisPop.saf.idx norNcen_HisPop.saf.idx  -fstout HisNorNCen_HisSouth -whichFst 1  -sfs HisNorNCen_HisSouth.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats HisNorNCen_HisSouth.fst.idx

