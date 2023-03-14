#Fst with ANGSD/v.921
################################################ Fst By Site ######################################################################################
#Number of bams in each population: 
# San Juan 5
# Coloso 3
# Cauca 1
# San Juan (Modern) 2
# Ciebal (Modern) 4
# Tierralta 5
# Planeta Rica 4
# Caracas 1
# Unknown 1
# Mutata 3
# Turbo 2
# Ajorna 1
# Tuplena (Modern) 2

#make a bash file per population with the following 
#run with  sbatch -c 5 --mem-per-cpu 4G --
# -P is number of threads
# 1. Estimate GLs

cd $STATS/Fst/Fst_site

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Arjona_HisPop.txt -doSaf 1 -out Arjona_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Caracas_HisPop.txt -doSaf 1 -out Caracas_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Cauca_HisPop.txt -doSaf 1 -out Cauca_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Ciebal_ModPop.txt  -doSaf 1 -out Ciebal_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Coloso_HisPop.txt -doSaf 1 -out Coloso_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Mutata_HisPop.txt -doSaf 1 -out Mutata_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam PlanetaRica_HisPop.txt -doSaf 1 -out PlanetaRica_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam SanJuan_HisPop.txt -doSaf 1 -out SanJuan_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam SanJuan_ModPop.txt -doSaf 1 -out SanJuan_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Tierralta_HisPop.txt -doSaf 1 -out Tierralta_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Turbo_HisPop.txt -doSaf 1 -out Turbo_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20

/projects/mjolnir1/apps/conda/angsd-0.921/bin/angsd -bam Tuplena_ModPop.txt -doSaf 1 -out Tuplena_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20


# 2. Estimate 2-dimensional frequency spectrum from the site allele frequency likelihoods 
#create bash files with: 
#run with sbatch -c 1 --mem-per-cpu 4G -J Fst -- 

# Historical Arjona vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Caracas_HisPop.saf.idx >Ajrona_Caracas.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Cauca_HisPop.saf.idx >Ajrona_Cauca.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Tierralta_HisPop.saf.idx >Ajrona_Tier.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Coloso_HisPop.saf.idx >Ajrona_Coloso.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Mutata_HisPop.saf.idx >Ajrona_Mutata.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx PlanetaRica_HisPop.saf.idx >Ajrona_PlantRica.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx SanJuan_HisPop.saf.idx >Ajrona_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Arjona_HisPop.saf.idx Turbo_HisPop.saf.idx >Ajrona_Turbo.ml

# Historical Caracas vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx Cauca_HisPop.saf.idx > Caracas_Cauca.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx Coloso_HisPop.saf.idx > Caracas_Coloso.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx Mutata_HisPop.saf.idx > Caracas_Mutata.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx PlanetaRica_HisPop.saf.idx > Caracas_PlantRica.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx SanJuan_HisPop.saf.idx > Caracas_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx Tierralta_HisPop.saf.idx > Caracas_Tier.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Caracas_HisPop.saf.idx Turbo_HisPop.saf.idx > Caracas_Turbo.ml

# Historical Cauca vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx Coloso_HisPop.saf.idx > Cauca_Coloso.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx Mutata_HisPop.saf.idx > Cauca_Mutata.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx PlanetaRica_HisPop.saf.idx > Cauca_PlantRica.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx SanJuan_HisPop.saf.idx > Cauca_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx Tierralta_HisPop.saf.idx > Cauca_Tier.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Cauca_HisPop.saf.idx Turbo_HisPop.saf.idx > Cauca_Turbo.ml

# Historical Coloso vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Coloso_HisPop.saf.idx Mutata_HisPop.saf.idx > Coloso_Mutata.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Coloso_HisPop.saf.idx PlanetaRica_HisPop.saf.idx > Coloso_PlantRica.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Coloso_HisPop.saf.idx SanJuan_HisPop.saf.idx > Coloso_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Coloso_HisPop.saf.idx Tierralta_HisPop.saf.idx > Coloso_Teir.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Coloso_HisPop.saf.idx Turbo_HisPop.saf.idx > Coloso_Turbo.ml

# Historical Mutata vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Mutata_HisPop.saf.idx PlanetaRica_HisPop.saf.idx > Mutata_PlantRica.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Mutata_HisPop.saf.idx SanJuan_HisPop.saf.idx > Mutata_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Mutata_HisPop.saf.idx Tierralta_HisPop.saf.idx > Mutata_Teir.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Mutata_HisPop.saf.idx Turbo_HisPop.saf.idx > Mutata_Turbo.ml

# Historical Planeta Rica vs all historical
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS PlanetaRica_HisPop.saf.idx SanJuan_HisPop.saf.idx > PlantRica_SJN.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS PlanetaRica_HisPop.saf.idx Tierralta_HisPop.saf.idx > PlantRica_Tier.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS PlanetaRica_HisPop.saf.idx Turbo_HisPop.saf.idx > PlantRica_Turbo.ml

# Historical San Juan vs all historical and Modern San Juan 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS SanJuan_HisPop.saf.idx Tierralta_HisPop.saf.idx > SJNHis_Teir.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS SanJuan_HisPop.saf.idx Turbo_HisPop.saf.idx > SJNHis_Turbo.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS SanJuan_HisPop.saf.idx SanJuan_ModPop.saf.idx > SJNHis_SJNMod.ml

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS SanJuan_HisPop.saf.idx Ciebal_ModPop.saf.idx > SJNHis_CEIMod.ml


# Historical Tierralta vs all historical 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Tierralta_HisPop.saf.idx Turbo_HisPop.saf.idx > Teir_Turbo.ml

# Historical Turbo vs Modrn Tuplena
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS Turbo_HisPop.saf.idx Tuplena_ModPop.saf.idx > Turbo_Tupelena.ml



# 3. pops fsts - send all below through sbtach --> sbatch -c 1 --mem-per-cpu 4G -J 
# the result file (*.Fst) contains two columns FST.Unweight(column 1) and Fst.Weight (Column 2)
#The weighted fst for a region is the ratio between the sum of As and the sum of B. The unweighted is the mean of the persite ratios.
# use the weighted Fst becuase it takes into account the the various sizes of the populations

## Historical:  Ajorna 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Caracas_HisPop.saf.idx -fstout Ajrona_Caracas -whichFst 1 -sfs Ajrona_Caracas.ml &
wait 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Caracas.fst.idx > Ajrona_Caracas.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Cauca_HisPop.saf.idx -fstout Ajrona_Cauca -whichFst 1  -sfs Ajrona_Cauca.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Cauca.fst.idx  > Ajrona_Cauca.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Ajrona_Tier -whichFst 1  -sfs Ajrona_Tier.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Tier.fst.idx > Ajrona_Tier.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Coloso_HisPop.saf.idx -fstout Ajrona_Coloso -whichFst 1  -sfs Ajrona_Coloso.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Coloso.fst.idx > Ajrona_Coloso.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Mutata_HisPop.saf.idx -fstout Ajrona_Mutata -whichFst 1  -sfs Ajrona_Mutata.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Mutata.fst.idx > Ajrona_Mutata.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx PlanetaRica_HisPop.saf.idx -fstout Ajrona_PlantRica -whichFst 1  -sfs Ajrona_PlantRica.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_PlantRica.fst.idx > Ajrona_PlantRica.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout Ajrona_SJN -whichFst 1  -sfs Ajrona_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_SJN.fst.idx > Ajrona_SJN.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Arjona_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout Ajrona_Turbo -whichFst 1  -sfs Ajrona_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Ajrona_Turbo.fst.idx > Ajrona_Turbo.Fst


# Historical Caracas vs all historical

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx Cauca_HisPop.saf.idx -fstout Caracas_Cauca -whichFst 1  -sfs Caracas_Cauca.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_Cauca.fst.idx > Caracas_Cauca.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx Coloso_HisPop.saf.idx -fstout Caracas_Coloso -whichFst 1  -sfs Caracas_Coloso.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_Coloso.fst.idx > Caracas_Coloso.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx Mutata_HisPop.saf.idx -fstout Caracas_Mutata -whichFst 1  -sfs Caracas_Mutata.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_Mutata.fst.idx > Caracas_Mutata.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx PlanetaRica_HisPop.saf.idx -fstout Caracas_PlantRica -whichFst 1  -sfs Caracas_PlantRica.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_PlantRica.fst.idx > Caracas_PlantRica.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout Caracas_SJN -whichFst 1  -sfs Caracas_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_SJN.fst.idx > Caracas_SJN.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Caracas_Tier -whichFst 1  -sfs Caracas_Tier.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_Tier.fst.idx > Caracas_Tier.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Caracas_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout Caracas_Turbo -whichFst 1  -sfs Caracas_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Caracas_Turbo.fst.idx > Caracas_Turbo.Fst

# Historical Cauca vs all historical

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx Coloso_HisPop.saf.idx -fstout Cauca_Coloso -whichFst 1  -sfs Cauca_Coloso.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_Coloso.fst.idx > Cauca_Coloso.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx Mutata_HisPop.saf.idx -fstout Cauca_Mutata -whichFst 1  -sfs Cauca_Mutata.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_Mutata.fst.idx > Cauca_Mutata.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx PlanetaRica_HisPop.saf.idx -fstout Cauca_PlantRica -whichFst 1  -sfs Cauca_PlantRica.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_PlantRica.fst.idx > Cauca_PlantRica.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout Cauca_SJN -whichFst 1  -sfs Cauca_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_SJN.fst.idx > Cauca_SJN.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Cauca_Tier -whichFst 1  -sfs Cauca_Tier.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_Tier.fst.idx > Cauca_Tier.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Cauca_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout Cauca_Turbo -whichFst 1  -sfs Cauca_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Cauca_Turbo.fst.idx > Cauca_Turbo.Fst

# Historical Coloso vs all historical

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Coloso_HisPop.saf.idx Mutata_HisPop.saf.idx -fstout Coloso_Mutata -whichFst 1  -sfs Coloso_Mutata.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Coloso_Mutata.fst.idx > Coloso_Mutata.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Coloso_HisPop.saf.idx PlanetaRica_HisPop.saf.idx -fstout Coloso_PlantRica -whichFst 1  -sfs Coloso_PlantRica.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Coloso_PlantRica.fst.idx > Coloso_PlantRica.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Coloso_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout Coloso_SJN -whichFst 1  -sfs Coloso_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Coloso_SJN.fst.idx > Coloso_SJN.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Coloso_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Coloso_Teir -whichFst 1  -sfs Coloso_Teir.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Coloso_Teir.fst.idx > Coloso_Teir.Fst

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Coloso_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout Coloso_Turbo -whichFst 1  -sfs Coloso_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Coloso_Turbo.fst.idx > Coloso_Turbo.Fst

# Historical Mutata vs all historical

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Mutata_HisPop.saf.idx PlanetaRica_HisPop.saf.idx -fstout Mutata_PlantRica -whichFst 1  -sfs Mutata_PlantRica.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Mutata_PlantRica.fst.idx > Mutata_PlantRica.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Mutata_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout Mutata_SJN -whichFst 1  -sfs Mutata_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Mutata_SJN.fst.idx > Mutata_SJN.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Mutata_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Mutata_Teir -whichFst 1  -sfs Mutata_Teir.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Mutata_Teir.fst.idx > Mutata_Teir.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Mutata_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout Mutata_Turbo -whichFst 1  -sfs Mutata_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Mutata_Turbo.fst.idx > Mutata_Turbo.Fst


# Historical Planeta Rica vs all historical

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index PlanetaRica_HisPop.saf.idx SanJuan_HisPop.saf.idx -fstout PlantRica_SJN -whichFst 1  -sfs PlantRica_SJN.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats PlantRica_SJN.fst.idx > PlantRica_SJN.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index PlanetaRica_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout PlantRica_Tier -whichFst 1  -sfs PlantRica_Tier.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats PlantRica_Tier.fst.idx > PlantRica_Tier.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index PlanetaRica_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout PlantRica_Turbo -whichFst 1  -sfs PlantRica_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats PlantRica_Turbo.fst.idx > PlantRica_Turbo.Fst

# Historical San Juan vs all historical and Modern Pops
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index SanJuan_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout SJNHis_Teir -whichFst 1  -sfs SJNHis_Teir.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats SJNHis_Teir.fst.idx > SJN_Teir.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index SanJuan_HisPop.saf.idx Turbo_HisPop.saf.idx -fstout SJNHis_Turbo -whichFst 1  -sfs SJNHis_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats SJNHis_Turbo.fst.idx > SJNHis_Turbo.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index SanJuan_HisPop.saf.idx SanJuan_ModPop.saf.idx -fstout SJNHis_SJNMod -whichFst 1  -sfs SJNHis_SJNMod.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats SJNHis_SJNMod.fst.idx > SJNHis_SJNMod.Fst


/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index SanJuan_HisPop.saf.idx Ciebal_ModPop.saf.idx -fstout SJNHis_CEIMod -whichFst 1  -sfs SJNHis_CEIMod.ml &
wait 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats SJNHis_CEIMod.fst.idx > SJNHis_CEIMod.Fst

# Historical Tierralta vs all historical 
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Turbo_HisPop.saf.idx Tierralta_HisPop.saf.idx -fstout Teir_Turbo -whichFst 1  -sfs Teir_Turbo.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Teir_Turbo.fst.idx > Teir_Turbo.Fst

# Historical Turbo vs Modrn Tuplena
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index Turbo_HisPop.saf.idx Tuplena_ModPop.saf.idx -fstout Turbo_Tupelena -whichFst 1  -sfs Turbo_Tupelena.ml &
wait
/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Turbo_Tupelena.fst.idx > Turbo_Tupelena.Fst

# Combine all the Fst into one file with file name, FST.Unweight and Fst.Weight 
# the higher the Fst the more differentiated they are 

for file in *Fst
 do 
    filename=$file
 echo "${filename%.*}" >> file1.list
 done

cat *Fst >> file2.list

paste file1.list file2.list  > AllFst.txt
sed 's/ /\t/g' AllFst.txt > AllPopsFst.txt

################################################ By Region Fst ######################################################################################

#Number of bams in each population:
# 9 north_HisPop.txt (San Juan, Coloso, Cauca)
# 6 north_ModPop.txt (San Juan, Ciebal)
# 11 central_HisPop.txt (Tierralta, Planeta Rica, Caracas, Unknown)
# 6 south_HisPop.txt (Mutata, Turbo, Ajorna)
# 2 south_ModPop.txt (Tuplena)

#make a bash file per population with the following
#run with  sbatch -c 5 --mem-per-cpu 4G --
# -P is number of threads
# 1. Estimate GLs

#cd $STATS/Fst/Fst_region

#/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam central_HisPop.txt  -doSaf 1 -out central_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 5

#/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam north_HisPop.txt  -doSaf 1 -out north_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 4

#/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam south_HisPop.txt -doSaf 1 -out south_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 3

#/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam north_ModPop.txt -doSaf 1 -out north_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 3

#/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam south_ModPop.txt -doSaf 1 -out south_ModPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 2


# 2. Estimate 2-dimensional frequency spectrum from the site allele frequency likelihoods (for each pair, 3 in this case)
#create bash files with:

# Historical North vs central
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS north_HisPop.saf.idx central_HisPop.saf.idx >HisNorth_HisCentral.ml

# Historical North vs South
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS north_HisPop.saf.idx south_HisPop.saf.idx >HisNorth_HisSouth.ml

# Historical Central vs South
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS central_HisPop.saf.idx south_HisPop.saf.idx >HisSouth_HisCentral.ml

# Historical North vs Modern North
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS north_HisPop.saf.idx north_ModPop.saf.idx >HisNorth_ModNorth.ml

# Historical South vs Modern South
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS south_HisPop.saf.idx south_ModPop.saf.idx >HisSouth_ModSouth.ml

# Historical Central vs Modern South
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS central_HisPop.saf.idx south_ModPop.saf.idx >HisCentral_ModSouth.ml

# Historical Central vs Modern North
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS central_HisPop.saf.idx north_ModPop.saf.idx >HisCentral_ModNorth.ml

# Modern South vs Modern North
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS south_ModPop.saf.idx north_ModPop.saf.idx >ModSouth_ModNorth.ml

# His NorNcentral vs HisSouth
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS norNcen_HisPop.saf.idx south_HisPop.saf.idx >HisNorNCen_HisSouth.ml


# run the above with sbatch
#sbatch -c 4 --mem-per-cpu 20G --time=7-00 -J HisCentral_ModNorth_SFS  -- scripts/runSFS_HisCentral_ModNorth.sh


# 3. pops fsts - send all below through sbtach --> sbatch -c 1 --mem-per-cpu 4G --time=7-00 -J

## Historical:  North | Central | South
#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index north_HisPop.saf.idx central_HisPop.saf.idx  south_HisPop.saf.idx -fstout Histocial_only -whichFst 1 -sfs HisNorth_HisCentral.ml -sfs HisNorth_HisSouth.ml -sfs HisSouth_HisCentral.ml

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Histocial_only.fst.idx

## Modern only

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index south_ModPop.saf.idx north_ModPop.saf.idx  -fstout Modern_only -whichFst 1  -sfs Modern_only.ml

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats Modern_only.fst.idx

# Only north

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index north_HisPop.saf.idx north_ModPop.saf.idx  -fstout North_only -whichFst 1  -sfs HisNorth_ModNorth.ml

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats North_only.fst.idx


#only South

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index south_HisPop.saf.idx south_ModPop.saf.idx  -fstout South_only -whichFst 1  -sfs HisSouth_ModSouth.ml

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats South_only.fst.idx


## 2Regions historical only (North+Central vs South)

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index south_HisPop.saf.idx norNcen_HisPop.saf.idx  -fstout HisNorNCen_HisSouth -whichFst 1  -sfs HisNorNCen_HisSouth.ml

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats HisNorNCen_HisSouth.fst.idx


# 4. print stats on 50K windows --> sbatch -c 1 --mem-per-cpu 4G -J

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats2 Histocial_only.fst.idx -win 50000 -step 5000 >Histocial_only_50Kwin

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats2 North_only.fst.idx -win 50000 -step 5000 >North_only_50Kwin

#/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats2 South_only.fst.idx -win 50000 -step 5000 >South_only_50Kwin
