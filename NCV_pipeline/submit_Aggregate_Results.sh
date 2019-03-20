#!/bin/bash -l

while getopts ":f:" opt; do
  case $opt in
    f) input=$OPTARG;;
  esac
done
echo "$input"

workDir=`pwd`
baseName=`basename "$input" .vcf.gz`

exec &> $baseName'_.motifBreakR.log'

module load gcc/7.2.0
module load mysql/5.0.83
module load gsl/2.1
module load R/3.4.3

# use - fileName_vcf.gz (update file name)
Rscript Aggregate_Results.R -f $input
