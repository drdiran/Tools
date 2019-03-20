#!/bin/bash

##############################################################
# Developer: Oyediran Akinrinade
# Arguments: vcf file in the current wkdir
# Split WGS vcf file by chr
###############################################################

while getopts ":a:" opt; do
  case $opt in
    a) vcf=$OPTARG;;
  esac
done
echo -e "$vcf"
vcfFile=$vcf

# create analysis and results directories
workDir=`pwd`
mkdir DeepSEA_Analysis
mkdir DeepSEA_Result
mkdir DeepSEA_Analysis/SCRIPTS
mkdir DeepSEA_Analysis/LOG
scriptDir=$workDir/DeepSEA_Analysis/SCRIPTS
log=$workDir/DeepSEA_Analysis/LOG

# Tools dir
snpEff="~/TOOLS/snpEff"

coor2fasta="~/TOOLS/DeepSEA/0.94b/0_coor2fasta.R"
fasta2input="~/TOOLS/DeepSEA/0.94b/1_fasta2input.py"
DeepSEA="~/TOOLS/DeepSEA/0.94b/2_DeepSEA.lua"
h5ToOutput="~/TOOLS/DeepSEA/0.94b/3_h5ToOutput.py"
evo="~/TOOLS/DeepSEA/0.94b/evoevalues.production.py"
deepsea_cpu="~/TOOLS/DeepSEA/0.94b/deepsea.cpu"
deepsea_res="~/TOOLS/DeepSEA/0.94b/resources"
logistic="~/TOOLS/DeepSEA/0.94b/logistic"
GerpScore="~/TOOLS/DeepSEA/0.94/resources/Gerp"
phastConsScore="~/TOOLS/DeepSEA/0.94/resources/phastCons"
phyloPScore="~/TOOLS/DeepSEA/0.94/resources/phyloP"

module load snpEff/4.3
module load java/1.8.0_91
module load vep/90
module load bcftools/1.4

# create txt files for overall results
touch DeepSEA_Result/allChrom_alt.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.alt | head -1 >> DeepSEA_Result/allChrom_alt.txt
touch DeepSEA_Result/allChrom_diff.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.diff | head -1 >> DeepSEA_Result/allChrom_diff.txt
touch DeepSEA_Result/allChrom_evalue.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.evalue | head -1 >> DeepSEA_Result/allChrom_evalue.txt
touch DeepSEA_Result/allChrom_funsig.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.funsig | head -1 >> DeepSEA_Result/allChrom_funsig.txt
touch DeepSEA_Result/allChrom_logfoldchange.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.logfoldchange | head -1 >> DeepSEA_Result/allChrom_logfoldchange.txt
touch DeepSEA_Result/allChrom_ref.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.ref | head -1 >> DeepSEA_Result/allChrom_ref.txt
touch DeepSEA_Result/allChrom_snpclass.txt
cat ~/TOOLS/DeepSEA/infile.vcf.out.snpclass | head -1 >> DeepSEA_Result/allChrom_snpclass.txt
touch DeepSEA_Result/allChrom_summary.txt
touch DeepSEA_Result/allChrom_wt1100_fasta_ref_evoall.txt
touch DeepSEA_Result/allChrom_wt1100_fasta_ref_evo_evalues.txt

ln -s $workDir/$vcfFile ./DeepSEA_Analysis/

#ln -s $vcfFile ./DeepSEA_Analysis/

cd ./DeepSEA_Analysis
#gzip -c -d $vcfFile > "${FilebaseName}.vcf"
java -jar -Xmx20g $snpEff/SnpSift.jar split $vcfFile
#rm -rf "${FilebaseName}.vcf"

ln -s $coor2fasta .
ln -s $fasta2input .
ln -s $DeepSEA .
ln -s $h5ToOutput .
ln -s $evo .
ln -s $deepsea_cpu .
ln -s $logistic .
mkdir resources
ln -s $deepsea_res/* ./resources
ln -s $GerpScore ./resources
ln -s $phastConsScore ./resources
ln -s $phyloPScore ./resources

for f in *_gnomADgenome.chr*.vcf;
do
varCount=`cat $f | grep -v "#" | wc -l`
if ((varCount > 0)); then
baseName=`basename "$f" .vcf`
mkdir $baseName
echo -e "
#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=20g
#PBS -l walltime=24:00:00
#PBS -N DeepSea
#PBS -d $workDir/DeepSEA_Analysis
#PBS -o $log/$baseName.out.logic
#PBS -e $log/$baseName.err.log

module load DeepSEA/0.94b
module load R/3.4.3

rundeepsea.py ./$f ./$baseName
sed 1d ./$baseName/infile.vcf.out.alt >> ../DeepSEA_Result/allChrom_alt.txt
sed 1d ./$baseName/infile.vcf.out.diff >> ../DeepSEA_Result/allChrom_diff.txt
sed 1d ./$baseName/infile.vcf.out.evalue >> ../DeepSEA_Result/allChrom_evalue.txt
sed 1d ./$baseName/infile.vcf.out.funsig >> ../DeepSEA_Result/allChrom_funsig.txt
sed 1d ./$baseName/infile.vcf.out.logfoldchange >> ../DeepSEA_Result/allChrom_logfoldchange.txt
sed 1d ./$baseName/infile.vcf.out.ref >> ../DeepSEA_Result/allChrom_ref.txt
sed 1d ./$baseName/infile.vcf.out.snpclass >> ../DeepSEA_Result/allChrom_snpclass.txt
cat ./$baseName/infile.vcf.out.summary >> ../DeepSEA_Result/allChrom_summary.txt
cat ./$baseName/infile.vcf.wt1100.fasta.ref.vcf.evoall >> ../DeepSEA_Result/allChrom_wt1100_fasta_ref_evoall.txt
cat ./$baseName/infile.vcf.wt1100.fasta.ref.vcf.evo.evalues >> ../DeepSEA_Result/allChrom_wt1100_fasta_ref_evo_evalues.txt
rm -rf $f

" > $scriptDir/$baseName'.sh'

chmod +x $scriptDir/$baseName'.sh'
echo "$baseName"
qsub $scriptDir/$baseName'.sh'

else
  echo "File $f is empty"
  echo "nothing done"
fi

done
