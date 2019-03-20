#!/bin/bash -l


while getopts ":f:" opt; do
  case $opt in
    f) input=$OPTARG;;
  esac
done
echo "$input"

workDir=`pwd`
baseName=`basename "$input" .vcf.gz`


echo "
#!/bin/bash -l
#PBS -l nodes=1:ppn=8
#PBS -l vmem=40g
#PBS -l walltime=23:00:00
#PBS -N Aggregate
#PBS -d $workDir
#PBS -o $workDir/LOG/$baseName'.Aggregate.out.logic'
#PBS -e $workDir/LOG/$baseName'.Aggregate.err.logic'

module load gcc/7.2.0
module load mysql/5.0.83
module load gsl/2.1
module load R/3.4.3

# use - fileName_vcf.gz (update file name)
Rscript Aggregate_Results.R -f $input

" > $workDir/SCRIPTS/$baseName'_Aggregate.sh'
chmod +x $workDir/SCRIPTS/$baseName'_Aggregate.sh'
echo "$baseName"
qsub $workDir/SCRIPTS/$baseName'_Aggregate.sh'
