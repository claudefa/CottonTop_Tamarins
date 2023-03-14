#!/bin/bash

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst index north_ModPop.saf.idx norNcen_HisPop.saf.idx  -fstout HisNorNCen_ModNorth -whichFst 1  -sfs HisNorNCen_ModNorth.ml &

wait

/projects/mjolnir1/apps/conda/angsd-0.921/bin/realSFS fst stats HisNorNCen_ModNorth.fst.idx 

