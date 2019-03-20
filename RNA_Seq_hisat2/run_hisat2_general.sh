#!/bin/bash

##############################################################
#This script will create scripts to run the hisat2 pipeline on samples in one directory
#it will then submit the scripts as a job to the cluster
#index information: pre-built index built using splice sites and snps with version 75 of grch37
###############################################################

# Ref files
INDEX=~/hisat2_index/built_index/grch37_snp_tran/genome_snp_tran  #location of index files for hisat                                                                              
GTF=~/hisat2_index/Homo_sapiens.GRCh37.75.gtf  #annotation file needed for stringtie  

#declare the indir for the input files and where to output the files
fastqDir=`pwd`

echo -e "$fastqDir/TRIMMED_FASTQ\n$fastqDir/fastqc\n$fastqDir/SCRIPTS\n$fastqDir/ALIGNMENT\n$fastqDir/LOG\n$fastqDir/ballgown" > dir_list.txt

while read dir; do
    if [[ -d $dir ]]; then
        echo "Dir exists"
    else
        mkdir $dir 
    fi  
done < dir_list.txt

rm -rf dir_list.txt

scriptDir=$fastqDir/SCRIPTS
OUTDIR=$fastqDir/ALIGNMENT #directory for output directories of each sample
LOG=$fastqDir/LOG
transcriptDir=$OUTDIR/ballgown #output directory for stringtie
ballgownDir=$OUTDIR/ballgown

for f in `ls $fastqDir/*'_R1.fastq.gz' | sed 's/_R1.fastq.gz//' `
do
baseName=$(basename ${f})

echo -e "
#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=45g
#PBS -l walltime=24:00:00
#PBS -N hisat2
#PBS -o $LOG/$baseName.out.logic
#PBS -e $LOG/$baseName.err.log

module load fastqc
module load trim_galore/0.4.4  # uses cutadapt/1.8 as prereq
module load hisat
module load samtools/1.3.1
module load stringtie

# check qc with fastqc
fastqc $fastqDir/$baseName'_R1.fastq.gz' $fastqDir/$baseName'_R2.fastq.gz' -q -f fastq -o $fastqDir/fastqc

# Read trimming and adapter removal/cutting
trim_galore -q 20 --illumina $fastqDir/$baseName'_R1.fastq.gz' $fastqDir/$baseName'_R2.fastq.gz' --paired --gzip --length 25 -o $fastqDir/TRIMMED_FASTQ

# check qc with trimmed reads and run hisat
if [ -f $fastqDir/TRIMMED_FASTQ/$baseName'_R1_val_1.fq.gz' ]; then
	echo "performing qc on trimmed reads"
fastqc $fastqDir/TRIMMED_FASTQ/$baseName'_R1_val_1.fq.gz' -2 $fastqDir/TRIMMED_FASTQ/$baseName'_R2_val_2.fq.gz' -q -f fastq -o $fastqDir/TRIMMED_FASTQ
        echo "performing hisat2 alignment"  
hisat2 --dta -p 8 -x $INDEX -1 $fastqDir/TRIMMED_FASTQ/$baseName'_R1_val_1.fq.gz' -2 $fastqDir/TRIMMED_FASTQ/$baseName'_R2_val_2.fq.gz' -S $OUTDIR/$baseName'.sam'
samtools sort -@ 8 -o $OUTDIR/$baseName'.bam' $OUTDIR/$baseName'.sam' 

else
	echo "not able to perform hisat2 alingment"
fi

## add -A as parameter to get gene abundance (input into R)

if [ -f $OUTDIR/$baseName'.bam' ]; then
	rm $OUTDIR/$baseName'.sam'
	mkdir $ballgownDir/$baseName
	echo "determining the transcript abundance and gene counts"
	stringtie $OUTDIR/$baseName'.bam' -G $GTF -l $baseName -o $ballgownDir/$baseName/$baseName'.gtf' -p 8 -e
else
	echo "not able to assemble transcripts"
fi

if [ -f $ballgownDir/$baseName/$baseName'.gtf' ]; then
        rm $OUTDIR/$baseName'.bam'
        echo "DONE"
else
        echo "could not remove bam file"
fi


" > $scriptDir/$baseName'.sh'

chmod +x $scriptDir/$baseName'.sh'
echo "$baseName"
qsub $scriptDir/$baseName'.sh'

done
