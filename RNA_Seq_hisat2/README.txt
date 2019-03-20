########################################################
##### README file for RNA-Seq data analysis pipeline####
##### Developer: Oyediran Akinrinade, PhD ##############
########################################################

The pipeline has two parts:
- 1. run_hisat2_general.sh
- 2. Prepare_DE_Analysis_File.sh (calls prepDE.py script)

## run_hisat2_general.sh
This script prepares/generates scripts for submitting job(s) for all samples and
submit each to the cluster. The script generated for each sample performs the following:
 - qc check with fastqc
 - adaptor trimming / removing with trim galore
 - qc on trimmed reads
 - runs hisat2 alignment
 - Estimate gene and transcript abundance with stringtie
 - checks if basName.gtf file exist. This file is needed for the second part.

 ## Prepare_DE_Analysis_File.sh
 This script calls the python script prepDE.py that aggregates the gene abundance
 info for all the samples. This generates gene count matrix file (gene_count_matrix.csv)
 and a transcript count matrix file (transcript_count_matrix.csv) needed for
 differential expression analysis with DeSeq2 or any tool of choice that accepts
 a count matrix.
