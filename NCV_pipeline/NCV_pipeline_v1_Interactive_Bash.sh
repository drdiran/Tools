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
dbsnp_150="/hpf/largeprojects/smital/references/homo_sapiens/hg19/dbSNP_150/GATK/All_20170403.vcf.gz"
plugins_dir="/hpf/largeprojects/smital/references/homo_sapiens/VEP/Plugins"
FATHMM_MKL_db="/hpf/largeprojects/smital/references/homo_sapiens/VEP/Plugins/fathmm-MKL_Current.tab.gz"
#vep87_dir="/hpf/largeprojects/smital/TOOLS/ensembl-tools-release-87/scripts/variant_effect_predictor"
#cache_dir="/hpf/tools/centos6/vep/90/genome"
cache_dir="/hpf/tools/centos6/vep/92/cache"
Human_Gnome_ref="/home/akinrina/Oyediran/WRKDIR/TOF/ucsc.hg19.fasta"
ref="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
fasta_refseq="/hpf/largeprojects/smital/references/homo_sapiens/VEP/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
dbNSFP2_9_2="/hpf/largeprojects/smital/references/homo_sapiens/hg19/dbNSFP_2.9.2/dbNSFP2.9.2.txt.gz"


#
Human_Ref="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
Human_Ref_bwa="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/bwa/human_g1k_v37_decoy_mod.fasta"
INDELS="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.indels.b37_mod.vcf"
MILLS="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/Mills_and_1000G_gold_standard.indels.b37_mod.vcf"
DBSNP="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/dbsnp_138.b37_mod.vcf"
HAPMAP="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/hapmap_3.3.b37_mod.vcf"
OMNI="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/1000G_omni2.5.b37_mod.vcf"
SNPS="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.snps.high_confidence.b37_mod.vcf"
dbsnp_147="/hpf/largeprojects/smital/references/homo_sapiens/hg19/dbSNP_147/All_20160601.vcf.gz"


#dbsnp_150="/hpf/largeprojects/smital/references/homo_sapiens/hg19/dbSNP_150/GATK/All_20170403.vcf.gz"
plugins_dir="/hpf/largeprojects/smital/references/homo_sapiens/VEP/Plugins"
#vep87_dir="/hpf/largeprojects/smital/TOOLS/ensembl-tools-release-87/scripts/variant_effect_predictor"
#cache_dir="/hpf/tools/centos6/vep/90/genome"
Human_Gnome_ref="/home/akinrina/Oyediran/WRKDIR/TOF/ucsc.hg19.fasta"
ref="/hpf/largeprojects/smital/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
#fasta_refseq="/hpf/largeprojects/smital/references/homo_sapiens/VEP/homo_sapiens/87_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
dbNSFP2_9_2="/hpf/largeprojects/smital/references/homo_sapiens/hg19/dbNSFP_2.9.2/dbNSFP2.9.2.txt.gz"
GNOMAD_Genome="/hpf/largeprojects/smital/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Exome="/hpf/largeprojects/smital/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Genome_decomposed="/hpf/largeprojects/smital/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"
GNOMAD_Exome_decomposed="/hpf/largeprojects/smital/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"


echo "
#!/bin/bash -l

exec &> $baseName'.log'

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

bash /hpf/largeprojects/smital/TOOLS/NCV_pipeline/NCV_firstSteps.sh -a $vcfFile

" > $scriptDir/$baseName'.sh'

chmod +x $scriptDir/$baseName'.sh'
echo "$baseName"
bash $scriptDir/$baseName'.sh' & 
