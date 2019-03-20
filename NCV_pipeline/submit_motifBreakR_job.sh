#!/bin/bash -l

exec &> $baseName'.log'

while getopts ":r:" opt; do
  case $opt in
    r) input=$OPTARG;;
  esac
done
echo "$input"

workDir=`pwd`
baseName=`basename "$input" .bed`

exec &> $baseName'motifBreakR.log'
## Please modify -d and -o to CWD
#PBS -l nodes=1:ppn=8
#PBS -l vmem=40g
#PBS -l walltime=24:00:00
#PBS -N motifBreakR
#PBS -j oe
#PBS -d $workDir
#PBS -o $workDir

module load gcc/7.2.0
module load mysql/5.0.83
module load gsl/2.1
module load R/3.4.3

# use - fileName_vep90_gnomADgenome_Regulatory_SNVs_MotifbreakR_Input.bed (update file name)
Rscript motifBreakR_pipeline.R -r $input
