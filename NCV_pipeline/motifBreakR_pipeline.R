#!/usr/bin/env Rscript

library(getopt)
library(motifbreakR)
#library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(GenomicRanges)
library(IRanges)

### FUNCTIONS ######################################################################################
usage <- function(){
  usage.value <-
    '
  motifBreakR_pipeline.R
  -help/-h (help)
  -motif_File/-r (Output from NCV_pipeline)
  ';
  cat(usage.value);
}

params <- matrix(
  c(
    'help', 'h', 0, 'logical',
    'motif_File', 'r', 1, 'character'
  ),
  byrow = TRUE,
  ncol = 4
);

opt <- getopt(params);

if(is.null(opt$motif_File)){
  print("no input file provided");
  usage();
  stop(options("show.error.messages" = FALSE));
}


# read in the promoter information of the CMP genes
funcFindOverlap=function(qry,sbj){
  k=findOverlaps(qry,sbj)
  match_hit=data.frame(names(qry)[queryHits(k)],names(sbj)[subjectHits(k)], stringsAsFactors=F)
  names(match_hit)= c('gene','sbjID')
  qrydf=as.data.frame(qry)
  sbjdf=as.data.frame(sbj)
  sbjdf=sbjdf[rownames(sbjdf) %in% match_hit$sbjID,]
  Strand=character(); GeneID=character();Transcription_Start=numeric()
  for (j in rownames(sbjdf)){
    ex=match_hit[match_hit$sbjID==j,"gene"]
    st=qrydf[rownames(qrydf)==ex,"STRAND_INFO"]
    ts=qrydf[rownames(qrydf)==ex,"TSS"]
    GeneID=append(GeneID,ex)
    Strand=append(Strand,st)
    Transcription_Start=append(Transcription_Start,ts)
  }
  sbjdf=data.frame(sbjdf,GeneID,Transcription_Start,Strand,stringsAsFactors=F)
  sbjdf=subset(sbjdf,select=c("GeneID","Transcription_Start","Strand"))
  return(sbjdf)
}

# Read the input file
Inputfile=opt$motif_File
snps.mb.frombedAll=snps.from.file(file = opt$motif_File,search.genome = BSgenome.Hsapiens.UCSC.hg19, format = "bed")
fname=unlist(strsplit(Inputfile,split="_Input"))[1]
OutputName=paste(fname,"MotifBreakR_Scores.RData",sep="_")


#snps.mb.frombedAll=snps.from.file(file = "ALL_TOF_Patients_CHD1_RegVariants_MotifbreakR_Input.bed",search.genome = BSgenome.Hsapiens.UCSC.hg19, format = "bed")

resultsAll_hoco <- motifbreakR(snpList = snps.mb.frombedAll, filterp = TRUE,
                               pwmList = hocomoco,
                               threshold = 1e-4,
                               method = "ic",
                               bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                               BPPARAM = BiocParallel::bpparam())

resultsAll_homer <- motifbreakR(snpList = snps.mb.frombedAll, filterp = TRUE,
                                pwmList = homer,
                                threshold = 1e-4,
                                method = "ic",
                                bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                BPPARAM = BiocParallel::bpparam())

resultsAll_encode <- motifbreakR(snpList = snps.mb.frombedAll, filterp = TRUE,
                                 pwmList = encodemotif,
                                 threshold = 1e-4,
                                 method = "ic",
                                 bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 BPPARAM = BiocParallel::bpparam())

# Subset to exclude weak effects
resultsAll_hoco=resultsAll_hoco[resultsAll_hoco$effect=="strong",]
resultsAll_homer=resultsAll_homer[resultsAll_homer$effect=="strong",]
resultsAll_encode=resultsAll_encode[resultsAll_encode$effect=="strong",]

save(resultsAll_hoco,resultsAll_homer,resultsAll_encode,snps.mb.frombedAll,file=OutputName)



