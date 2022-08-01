#!/bin/bash

bed=$1;
size=$2;
slide=$3;

if [ -z "$bed" ] || [ -z "$size" ] || [ -z "$slide" ];then printf "### This builds windows ###\n\nUSAGE:\n\nbash makewindows.sh Assembly_scaffolds_autosomes.bed 100000 50000 \n\n";exit;fi
/apps/BEDTools/2.26.0/bin/windowMaker -b $1 -w $2 -s $3 | /apps/BEDTools/2.26.0/bin/sortBed -i /dev/stdin | bgzip > $2.bed.gz;tabix -p bed $2.bed.gz
