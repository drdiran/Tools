#!/bin/bash

##############################################################
#This script will prepare csv files (gene & transcript levels) for DE analysis
#it has to be run from the ../ballgown/TranscriptCounts directory
###############################################################

module load python/2.7.15

# first prepare sample file for script prep.py
ls *.gtf | cat | awk -v OFS='\t' '{split($0,a,"_"); print a[1],$0}' > SampleList.txt

python2 prepDE.py -i SampleList.tx
