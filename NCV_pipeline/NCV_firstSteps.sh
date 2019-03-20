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


# workDir
workDir=`pwd`

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

# Split/decompose, block, MNP, CLUMPED etc
vt decompose $vcfFile | vt decompose_blocksub -a - -o $outFile'_decomposed.vcf'

## Remove ID & INFO field
# Reformat the vcf file
bcftools annotate -x ID,INFO -O v $outFile'_decomposed.vcf' -o $outFile'_decomposed_reformedvcf.vcf'
java -Xmx5g -jar $snpEff/SnpSift.jar annotate -id -noInfo $dbsnp_150 $outFile'_decomposed_reformedvcf.vcf' > $outFile'_decomposed_reformedvcf_dbsnp150.vcf'


## Annotate with vep/92
vep --assembly GRCh37 --everything --refseq --no_progress --dir_plugins $plugins_dir --fork 4 --cache --dir $cache_dir --no_stats --offline --fasta $fasta_refseq \
-i $outFile'_decomposed_reformedvcf_dbsnp150.vcf' -o $outFile'_vep92.vcf' \
-custom $GNOMAD_Genome_decomposed,gnomADgVep,vcf,exact,0,AF,AC,AN,Hom \
-custom $GNOMAD_Exome_decomposed,gnomADeVep,vcf,exact,0,AF,AC,AN,Hom \
--flag_pick --plugin dbNSFP,$dbNSFP2_9_2,GERP++_RS,CADD_raw,CADD_phred,ExAC_AC \
--plugin LoF --vcf \
--plugin FATHMM_MKL,$FATHMM_MKL_db

java -Xmx10G -jar $snpEff/SnpSift.jar annotate -v -noId -info AF,AC,Hom,AN -name gnomADg_ $GNOMAD_Genome_decomposed $outFile'_vep92.vcf' > $outFile'_vep92_gnomADgenome.vcf'

# remove inline Files
rm -rf $outFile'_decomposed.vcf'
rm -rf $outFile'_decomposed_reformedvcf.vcf'
rm -rf $outFile'_decomposed_reformedvcf_dbsnp150.vcf'
rm -rf $outFile'_vep92.vcf'

# Extract variants overlapping know RegulatoryFeature
filter_vep -i $outFile'_vep92_gnomADgenome.vcf' \
--format vcf -o $outFile'_vep92_gnomADgenome_Regulatory.vcf' \
--filter "Consequence matches regulatory_region_variant" --only_matched

# Extract SNVs overlapping know RegulatoryFeature
filter_vep -i $outFile'_vep92_gnomADgenome.vcf' \
--format vcf -o $outFile'_vep92_gnomADgenome_Regulatory_SNVs.vcf' \
--filter "Consequence matches regulatory_region_variant" --filter "VARIANT_CLASS matches SNV" --only_matched

## Convert to flat file
# Combined
bcftools query -l $outFile'_vep92_gnomADgenome_Regulatory.vcf' | tr '\n' '\t' | awk '{print "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tCSQ\t" $0}' > $outFile'_vep92_gnomADgenome_Regulatory.txt'
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/CSQ[\t%GT]\n' $outFile'_vep92_gnomADgenome_Regulatory.vcf' >> $outFile'_vep92_gnomADgenome_Regulatory.txt'

# SNVs
bcftools query -l $outFile'_vep92_gnomADgenome_Regulatory_SNVs.vcf' | tr '\n' '\t' | awk '{print "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tCSQ\t" $0}' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs.txt'
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/CSQ[\t%GT]\n' $outFile'_vep92_gnomADgenome_Regulatory_SNVs.vcf' >> $outFile'_vep92_gnomADgenome_Regulatory_SNVs.txt'

# Generate input files for Regulomedb query and coverage analysis
##awk -v OFS='\t' 'NR>1{print $1,$2-1,$2}' $outFile'_vep92_gnomADgenome_Regulatory_SNVs.txt' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input.txt'
##awk -v OFS='\t' '{gsub(/^chr/,""); print}' $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input.txt' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input_nochr.txt'

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/CSQ[\t%GT]\n' $outFile'_vep92_gnomADgenome_Regulatory_SNVs.vcf' | awk -v OFS='\t' '{print $1,$2-1,$2}' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input.txt'
cat $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input.txt' | awk -v OFS='\t' '{gsub(/^chr/,""); print}' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input_nochr.txt'

# generate motifbreakR input files
awk -v OFS='\t' 'NR>1{a=$2-1;b=0; c="+"; print $1,a,$2,$1":"a":"$4":"$5,b,c}' $outFile'_vep92_gnomADgenome_Regulatory_SNVs.txt' > $outFile'_vep92_gnomADgenome_Regulatory_SNVs_MotifbreakR_Input.bed'

# Generate coverage info
## Coverage
module unload tabix/0.2.6
module load htslib/1.4.1
touch $outFile'_vep92_gnomADgenome_Regulatory_SNVs_Coverage.txt'
covFileList=~/references/homo_sapiens/hg19/gnomAD/genome_coverage/genome/Coverage_processed/Genome_Coverage_Files.list
input=$outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input_nochr.txt'
output=$outFile'_vep92_gnomADgenome_Regulatory_SNVs_Coverage.txt'
while read p; do
 tabix -p bed ~/references/homo_sapiens/hg19/gnomAD/genome_coverage/genome/Coverage_processed/$p -R $input >> $output
done <$covFileList

## Generate motifBreakR prediction
Rscript motifBreakR_pipeline.R -r $outFile'_vep92_gnomADgenome_Regulatory_SNVs_MotifbreakR_Input.bed'

 ## Generate RegulomeDB score
Rscript RegulomeDB_Query_Pipeline.R -r $outFile'_vep92_gnomADgenome_Regulatory_SNVs_RegDB_Input.txt'
