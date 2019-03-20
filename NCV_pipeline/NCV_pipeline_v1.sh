#!/bin/bash

while getopts ":a:" opt; do
  case $opt in
    a) vcf=$OPTARG;;
  esac
done
echo "$vcf"

# vcf file
vcfFile=$vcf
baseName=`basename "$vcfFile" .vcf.gz`
outFile=`basename "$vcfFile" .vcf.gz`

#baseName=`basename "$vcfFile" .vcf`

# workDir
workDir=`pwd`
mkdir SCRIPTS
mkdir LOG
scriptDir=$workDir/SCRIPTS
log=$workDir/LOG

# Reference Files
dbsnp_150="~/references/homo_sapiens/hg19/dbSNP_150/GATK/All_20170403.vcf.gz"
plugins_dir="~/references/homo_sapiens/VEP/Plugins"
FATHMM_MKL_db="~/references/homo_sapiens/VEP/Plugins/fathmm-MKL_Current.tab.gz"
cache_dir="~/references/vep/92/cache"
Human_Gnome_ref="~/references/ucsc.hg19.fasta"
ref="~/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
fasta_refseq="~/references/homo_sapiens/VEP/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
dbNSFP2_9_2="~/references/homo_sapiens/hg19/dbNSFP_2.9.2/dbNSFP2.9.2.txt.gz"
Human_Ref="~/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
Human_Ref_bwa="~/references/homo_sapiens/v37_decoy/bwa/human_g1k_v37_decoy_mod.fasta"
INDELS="~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.indels.b37_mod.vcf"
MILLS="~/references/homo_sapiens/v37_decoy/gatk/Mills_and_1000G_gold_standard.indels.b37_mod.vcf"
DBSNP="~/references/homo_sapiens/v37_decoy/gatk/dbsnp_138.b37_mod.vcf"
HAPMAP="~/references/homo_sapiens/v37_decoy/gatk/hapmap_3.3.b37_mod.vcf"
OMNI="~/references/homo_sapiens/v37_decoy/gatk/1000G_omni2.5.b37_mod.vcf"
SNPS="~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.snps.high_confidence.b37_mod.vcf"
dbsnp_147="~/references/homo_sapiens/hg19/dbSNP_147/All_20160601.vcf.gz"
GNOMAD_Genome="~/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Exome="~/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Genome_decomposed="~/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"
GNOMAD_Exome_decomposed="~/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"

# Tools dir
snpEff="~/TOOLS/snpEff/4.3"

echo "
#!/bin/bash -l
#PBS -l nodes=1:ppn=8
#PBS -l vmem=40g
#PBS -l walltime=24:00:00
#PBS -N Vep_92
#PBS -d $workDir
#PBS -o $log/$baseName.out.logic
#PBS -e $log/$baseName.err.log

module unload tabix/0.2.6
module load DeepSEA/0.94b
module load htslib/1.4.1
module load vt
module load vcftools
module load bcftools
module load snpEff/4.3
module load vep/92
module load bedtools
module load gcc/7.2.0
module load mysql/5.0.83
module load gsl/2.1
module load R/3.4.3

# vcf file
vcfFile=$vcf
baseName=`basename "$vcfFile" .vcf.gz`
outFile=`basename "$vcfFile" .vcf.gz`

bash NCV_pipeline/NCV_firstSteps.sh -a $vcfFile

" > $scriptDir/$baseName'.sh'

chmod +x $scriptDir/$baseName'.sh'
echo "$baseName"
qsub $scriptDir/$baseName'.sh'
