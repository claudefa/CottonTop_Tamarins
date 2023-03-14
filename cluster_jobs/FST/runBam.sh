#!/bin/bash
NUCLEARGENOME=/projects/mjolnir1/people/vbz170/projects/CTT/Ref_Genome/Saguinus_midas_Full_Genome/NCBI_version/SaguinusMidas_NCBI.fasta

/projects/mjolnir1/apps/conda/angsd-0.935/bin/angsd -bam norNcen_HisPop.txt  -doSaf 1 -out norNcen_HisPop -anc $NUCLEARGENOME -gl 2 -P 5 -minMapQ 20 -minQ 20 -minInd 5

