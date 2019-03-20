#!/usr/bin/env Rscript
library(getopt)

### FUNCTIONS ######################################################################################
usage <- function(){
  usage.value <-
    '
  Aggregate_Results.R
  -help/-h (help)
  -fileName/-f (Starting vcf file)
  ';
  cat(usage.value);
}

params <- matrix(
  c(
    'help', 'h', 0, 'logical',
    'fileName', 'f', 1, 'character'
  ),
  byrow = TRUE,
  ncol = 4
);

opt <- getopt(params);

if(is.null(opt$fileName)){
  print("no input file provided");
  usage();
  stop(options("show.error.messages" = FALSE));
}

source("./Non_coding_variants_Analysis_Functions.R")
vep92Header="~/Non_Coding_Variants_Analysis_Resources/Vep92Annotation_header.txt"

# check if required files are present
baseName=unlist(strsplit(opt$fileName,split=".vcf"))[1]
# regulatory SNV file
RegSNV=paste(baseName,"_vep92_gnomADgenome_Regulatory_SNVs.txt",sep="")
# Coverage file
RegCov=paste(baseName,"_vep92_gnomADgenome_Regulatory_SNVs_Coverage.txt",sep="")
# RegulomeDB score
RegDBscore=paste(baseName,"_vep92_gnomADgenome_Regulatory_SNVs_RegDB_RegDB_Scores.RData",sep="")
#MotifBreakR Score
motifBreakRscore=paste(baseName,"_vep92_gnomADgenome_Regulatory_SNVs_MotifbreakR_MotifBreakR_Scores.RData",sep="")
# DeepSEA score
deepSea="./DeepSEA_Result/allChrom_funsig.txt"

reqFiles=c(RegSNV,RegCov,RegDBscore,motifBreakRscore,deepSea)

if(is.element("FALSE",file.exists(reqFiles))==T){
  print("Required file(s) not found");
  stop(options("show.error.messages" = FALSE));
} else {
  # Read in processed vcf file converted to txt file for SNVs
  allDat=read.delim(RegSNV,header=T,sep="\t",stringsAsFactors = F)
  # remove the last column
  allDat=allDat[,-ncol(allDat)]
  vepHeader=read.delim(vep92Header,header=F, sep="\t",stringsAsFactors = F)
  allDatVar=allDat[,1:8]
  genotype=allDat[,-c(1:8)]
  if (ncol(genotype)>1){
    Samples=do.call(rbind,lapply(1:nrow(genotype),FuncSum,y=genotype))
    allDatVar=cbind(allDatVar,Samples)
    GenoSummary=do.call(rbind,lapply(1:nrow(genotype),funcSummerizeGenotype_NC,PtGenotype=genotype))
    allDatVar=cbind(allDatVar,GenoSummary)
  } else {
    allDatVar$Zygosity=do.call(rbind,lapply(genotype,funcZygosity))
  }
  #split CSQ field
  CSQsplit=do.call(rbind,lapply(allDatVar$CSQ,funcSPlitCSQ))
  colnames(CSQsplit)=vepHeader$V1
  CSQsplit_2=do.call(rbind,lapply(allDatVar$CSQ,funcSPlitCSQextractSecondFeature))
  colnames(CSQsplit_2)=c("Feature_2","BIOTYPE_2")
  CSQsplitNeed=cbind(CSQsplit,CSQsplit_2)
  allDatVar=cbind(allDatVar,CSQsplitNeed)

  # Update varID
  varID=paste(allDatVar$REF,allDatVar$ALT,sep="/")
  allDatVar$varID=paste(allDatVar$CHROM,allDatVar$POS,sep=":")
  allDatVar$varID2=paste(allDatVar$CHROM,allDatVar$POS,varID,sep=":")

  # Exclude duplicates
  allDatVar=allDatVar[!duplicated(allDatVar$varID2),]

  # Update nearby TSS and gene info; Regulatory Feature Info
  sbj=allDatVar[,c("CHROM","POS","REF","ALT","varID")];sbj$varID=paste(sbj$varID,sbj$REF,sep=":"); sbj$varID=paste(sbj$varID,sbj$ALT,sep="/");
  rownames(sbj)=sbj$varID; sbj$end=sbj$POS; sbj=sbj[,c("CHROM","POS","end","varID")];colnames(sbj)=c("chr","start","end","Variant")
  sbj=with(sbj,GRanges(chr,IRanges(start,end,names=Variant),Variant=Variant,varPos=start))

  qry_enh=EnhancerUniverse_Compedium
  qry_pro=PromoterUniverse_Compedium
  allDatVar_Annot_enh=funcFindOverlap(qry_enh,sbj)
  allDatVar_Annot_pro=funcFindOverlap(qry_pro,sbj)

  colnames(allDatVar_Annot_enh)[-1]=paste(colnames(allDatVar_Annot_enh)[-1],"Enhancer",sep="_");
  colnames(allDatVar_Annot_pro)[-1]=paste(colnames(allDatVar_Annot_pro)[-1],"Promoter",sep="_");

  # Merge enhancer and promoter annotations with allDatVar
  allDatVar=merge(allDatVar,allDatVar_Annot_enh,by="varID2",all.x=T,sort=F)
  allDatVar=merge(allDatVar,allDatVar_Annot_pro,by="varID2",all.x=T,sort=F)

  # Update Gene ID
  allDatVar$Gene=allDatVar$NearbyGene_Enhancer
  allDatVar$Gene[which(is.na(allDatVar$NearbyGene_Enhancer)==T)]=allDatVar$NearbyGene_Promoter[which(is.na(allDatVar$NearbyGene_Enhancer)==T)]


  # Attach cov info from Gnomad
  cov=read.delim(RegCov,header=F,sep="\t", stringsAsFactors = F)
  cov$V1=paste("chr",cov$V1,sep="")
  cov$varID=paste(cov$V1,cov$V3,sep=":"); cov=cov[,c(6,4,5)]; colnames(cov)[c(2,3)]=c("meanCov","medianCov")
  cov=cov[!duplicated(cov$varID),]

  allDatVar=merge(allDatVar,cov,by="varID",sort=F,all.x=T)


  # load RegDB score - RData from R query
  load(RegDBscore,verbose=T)
  RegScore$V1=paste("chr",RegScore$V1,sep="")
  RegScore$varID=paste(RegScore$V1,RegScore$V3,sep=":")

  RegScoreAll=RegScore[!duplicated(RegScore$varID),]
  RegScoreAll=RegScoreAll[,c(6,5)]

  # append regulomedb score
  allDatVar=merge(allDatVar,RegScoreAll,by="varID",sort=F,all.x=T)

  # Add MotifBreakR predictions
  load(motifBreakRscore,verbose=T)
  # use only encode motifs
  resultsAll_encode$MotifEffect=as.character(do.call(rbind,lapply(1:length(resultsAll_encode),funcMotifEffect,y=resultsAll_encode)))

  resultsAll_encodedf=as.data.frame(resultsAll_encode,row.names = NULL)
  resultsAll_encodedf$varID2=paste(resultsAll_encodedf$seqnames,resultsAll_encodedf$snpPos,resultsAll_encodedf$REF,sep=":")
  resultsAll_encodedf$varID2=paste(resultsAll_encodedf$varID2,resultsAll_encodedf$ALT,sep="/")

  # Reduce resultsAll_encodedf to unique variant
  resultsAll_encodedfReduced=do.call(rbind,lapply(unique(resultsAll_encodedf$varID2),funcReduceMBR,y=resultsAll_encodedf))

  # update cardiacTF
  resultsAll_encodedfReduced$CardiacTF=do.call(rbind,lapply(1:nrow(resultsAll_encodedfReduced),funcUpdateCardiacTF,y=resultsAll_encodedfReduced,z=cardiacTF))

  allDatVar=merge(allDatVar,resultsAll_encodedfReduced,by="varID2",all.x=T,sort=F)

  # Update DeepSEA prediction score
  ## Read in DeepSEA fxnal signif score
  deepsea=read.delim(deepSea,header=T,sep=",",stringsAsFactors = F)
  deepsea=deepsea[,-1]
  deepsea$varID2=paste(deepsea$chr,deepsea$pos,deepsea$ref,sep=":")
  deepsea$varID2=paste(deepsea$varID2,deepsea$alt,sep="/")
  deepseaMap=deepsea[,c("varID2","Functional.significance.score")]; colnames(deepseaMap)[2]="DeepSEA_pvalue"

  #append deepsea score to allDatVar
  allDatVar=merge(allDatVar,deepseaMap,by="varID2",all.x=T,sort=F)

  # Add regulatory feature for LV
  #load("Ensembl_RegulatoryFeature_HeartActivity.RData",verbose = T)
  allDatVar=merge(allDatVar,Heart_Activity,by="Feature",all.x=T,sort=F)

  # Update
  allDatVar$varID=paste(allDatVar$varID,allDatVar$REF,sep="_"); allDatVar$varID=paste(allDatVar$varID,allDatVar$ALT,sep="/")

  # Summarize NCV scores
  Tier_score=do.call(rbind,lapply(1:nrow(allDatVar),funcReturnValue,q=allDatVar))
  allDatVar=cbind(allDatVar,Tier_score)
  outName=paste(baseName,"_Annotated",sep="")
  write.table(allDatVar,file=paste(outName,Sys.Date(),"Plain.txt",sep="_"),row.names=F,quote=F,sep="\t")
  save(allDatVar,file=paste(outName,Sys.Date(),"RFile.RData",sep="_"))
}
