#!/usr/bin/env Rscript

library(GenomicRanges)
library(IRanges)

load("~/Non_Coding_Variants_Analysis_Resources/Compedium_Enhancer_Promoter.RData",verbose = T)
load("~/Non_Coding_Variants_Analysis_Resources/CMPgenes_Enhancer_Promoter_Map_Compedium.RData",verbose=T)
load("~/Non_Coding_Variants_Analysis_Resources/Ensembl_RegulatoryFeature_HeartActivity.RData",verbose = T)
load("~/Non_Coding_Variants_Analysis_Resources/CHD_Genes_Enhancer_Promoter_Map.RData",verbose = T)
load("~/Non_Coding_Variants_Analysis_Resources/Overall_Cardiac_Regulome_Map.RData",verbose = T)

cardiacTF="~/Non_Coding_Variants_Analysis_Resources/Cardiac_TFs.txt"
vep92Header="~/Non_Coding_Variants_Analysis_Resources/Vep92Annotation_header.txt"


FuncSum=function(x,y){
  idt=y[x,];idt=data.frame(t(idt),stringsAsFactors = F);colnames(idt)="GT"
  idt$GAT=array("RAT",nrow(idt))
  RN=rownames(idt[grepl("^\\./.|0/0",idt$GT)==F,]);noRN=length(RN)
  RN=paste(RN,sep=", ",collapse=", ");mat=matrix(NA,0,2);colnames(mat)=c("SampleID","noSamples")
  mat=rbind(mat,c(RN,noRN));mat=as.data.frame(mat,stringsAsFactors=F)
  return(mat)
}


###
funcCollapseNames=function(y){
  if (length(y)==0){
    a=""
  } else {
    a=paste(y,sep=", ",collapse=", ")
  }
  return(a)
}

funcSummerizeGenotype_NC=function(z,PtGenotype){
  x=PtGenotype
  x=x[z,]
  y=x
  yt=data.frame(t(y),stringsAsFactors=F);colnames(yt)="genoType"
  yt$geno=do.call(rbind,strsplit(yt$genoType,split="\\:"))[,1]
  hom_var=yt[yt$geno %in%c("1/1","2/2","3/3"),]
  het_var=yt[yt$geno %in%c("0/1","1/0","0/2","2/0","0/3","3/0"),]
  heterozygotes=rownames(het_var);het_counts=length(heterozygotes)
  homozygotes=rownames(hom_var);hom_counts=length(homozygotes)
  heterozygotes=funcCollapseNames(heterozygotes)
  homozygotes=funcCollapseNames(homozygotes)
  df=data.frame(heterozygotes,het_counts,homozygotes,hom_counts,stringsAsFactors=F)
  return(df)
}

# zygosity function for a single sample input file
funcZygosity=function(x){
  if(x=="1/1" | x=="2/2" | x=="3/3"){
    Zygosity="Hom"
  } else {
    Zygosity="Het"
  }
  return(Zygosity)
}

## For 1kg data
FuncSum_1kg=function(x,y){
  idt=y[x,];idt=data.frame(t(idt),stringsAsFactors = F);colnames(idt)="GT"
  idt$GAT=array("RAT",nrow(idt))
  RN=rownames(idt[grepl("\\.\\|\\.|0\\|0",idt$GT)==F,]);noRN=length(RN)
  RN=paste(RN,sep=", ",collapse=", ");mat=matrix(NA,0,2);colnames(mat)=c("SampleID","noSamples")
  mat=rbind(mat,c(RN,noRN));mat=as.data.frame(mat,stringsAsFactors=F)
  return(mat)
}

###
funcCollapseNames_1kg=function(y){
  if (length(y)==0){
    a=""
  } else {
    a=paste(y,sep=", ",collapse=", ")
  }
  return(a)
}

funcSummerizeGenotype_1kg=function(z,PtGenotype){
  x=PtGenotype
  x=x[z,]
  y=x #;ynames=do.call(rbind,strsplit(colnames(y),split="\\."))[,4]
  #ynames=gsub("_","-",ynames); colnames(y)=ynames
  yt=data.frame(t(y),stringsAsFactors=F);colnames(yt)="genoType"
  yt$geno=do.call(rbind,strsplit(yt$genoType,split="\\:"))[,1]
  #ref=x$REF;alt=x$ALT
  hom_var=yt[yt$geno %in%c("1|1","2|2"),]
  het_var=yt[yt$geno %in%c("0|1","1|0","0|2","2|0"),]
  heterozygotes=rownames(het_var);het_counts=length(heterozygotes)
  homozygotes=rownames(hom_var);hom_counts=length(homozygotes)
  heterozygotes=funcCollapseNames_1kg(heterozygotes)
  homozygotes=funcCollapseNames_1kg(homozygotes)
  df=data.frame(heterozygotes,het_counts,homozygotes,hom_counts,stringsAsFactors=F)
  return(df)
}

###


funcUpdateVT=function(x){
  if(grepl("SNV",x)==T){
    a="SNV"
  } else {
    a="INDEL"
  }
  return(a)
}

funcUpdateRegFeat=function(x){
  if(grepl("Regulatory|regulatory",x)==T){
    a="Regulatory"
  } else {
    a=""
  }
  return(a)
}


# split CSQ/INFO
funcSPlitCSQ_vep90=function(x){
  mat=matrix(NA,0,87)
  if(grepl(",",x)==T){
    x=unlist(strsplit(x,split=","))[1]
  }
  xdt=unlist(strsplit(x,split="\\|"))
  if(length(xdt)==86){
    xdt=append(xdt,"")
  }
  mat=rbind(mat,xdt); mat=data.frame(mat,stringsAsFactors = F);rownames(mat)=NULL
  return(mat)
}

## split CSQ/INFO - vep92
funcSPlitCSQ=function(x){
  mat=matrix(NA,0,91)
  if(grepl(",",x)==T){
    x=unlist(strsplit(x,split=","))[1]
  }
  xdt=unlist(strsplit(x,split="\\|"))
  if(length(xdt)==90){
    xdt=append(xdt,"")
  }
  mat=rbind(mat,xdt); mat=data.frame(mat,stringsAsFactors = F);rownames(mat)=NULL
  return(mat)
}



funcSPlitCSQextractSecondFeature_1=function(x){
  mat=matrix(NA,0,2)
  if(grepl(",",x)==T){
    x=unlist(strsplit(x,split=","))[2]
    xdt=unlist(strsplit(x,split="\\|"))
    if(length(xdt)==90){
      xdt=append(xdt,"")
    }
    xdt=xdt[7:8]
  } else {
    xdt=c("","")
  }
  mat=rbind(mat,xdt); mat=data.frame(mat,stringsAsFactors = F);rownames(mat)=NULL
  return(mat)
}

funcSPlitCSQextractSecondFeature=function(x){
  mat=matrix(NA,0,2)
  if(grepl(",",x)==T){
    x=unlist(strsplit(x,split=","))[2]
    xdt=unlist(strsplit(x,split="\\|"))
    xdt=xdt[7:8]
  } else {
    xdt=c("","")
  }
  mat=rbind(mat,xdt); mat=data.frame(mat,stringsAsFactors = F);rownames(mat)=NULL
  return(mat)
}




# Function that generates nearby gene, TSS, distance, FeatureID
#qry=EnhancerUniverse # switch
#qry=EnhancerUniverse_Compedium
#qry=Promoter_Merged_Universe

#funcUP_DownStream=function(x,y){
#  strnd=y$Strand[x]; dss=y$Distance_fromTSS[x]
#  if(strnd=="+" & dss>0){
#    y="Downstream"
#  } else if (strnd=="+" & dss<0){
#    y="Upstream"
#  } else if (strnd=="-" & dss>0){
#    y="Upstream"
#  } else if (strnd=="-" & dss<0){
#    y="Downstream"
#  } else{
#    y="TSS"
#  }
#  return(y)
#}

funcUP_DownStream=function(x,y){
  strnd=y$Strand[x]; dss=y$Distance_fromTSS[x]
  if(is.na(strnd)==T | is.na(dss)==T) {
    y=""
  } else if(strnd=="+" & dss>0){
    y="Downstream"
  } else if (strnd=="+" & dss<0){
    y="Upstream"
  } else if (strnd=="-" & dss>0){
    y="Upstream"
  } else if (strnd=="-" & dss<0){
    y="Downstream"
  } else{
    y="TSS"
  }
  return(y)
}

funcFindOverlap=function(qry,sbj){
  k=findOverlaps(qry,sbj)
  match_hit=data.frame(names(qry)[queryHits(k)],names(sbj)[subjectHits(k)], stringsAsFactors=F)
  names(match_hit)= c('FeatureID','Variant')
  qrydf=as.data.frame(qry); rownames(qrydf)=NULL
  qrydf=qrydf[,c("GeneID","TSS","FeatureID","Strand")]
  match_hit=merge(match_hit,qrydf,by="FeatureID",all.x=T,sort=F)
  sbjdf=as.data.frame(sbj); sbjdf$SN=seq(1,nrow(sbjdf))
  sbjdf=merge(sbjdf,match_hit,by="Variant",all.x=T,sort=F)
  sbjdf$Distance_fromTSS=sbjdf$varPos-sbjdf$TSS
  sbjdf=sbjdf[order(sbjdf$SN,decreasing = F),]
  sbjdf$Strand=gsub("-1","-",sbjdf$Strand); sbjdf$Strand=gsub("1","+",sbjdf$Strand)
  sbjdf$varLocation=do.call(rbind,lapply(1:nrow(sbjdf),funcUP_DownStream,y=sbjdf))
  sbjdf=sbjdf[,c("Variant","GeneID","FeatureID","TSS","Strand","Distance_fromTSS","varLocation")]
  colnames(sbjdf)=c("varID2","NearbyGene","GeneFeatureID","NearbyTSS","geneStrand","Distance_fromTSS","varLocation")
  return(sbjdf)
}

## function for motifBreakR output: requires either promoter/enhancer region defn including GeneID and TSS
# Promoter_Merged_Universe or EnhancerUniverse as the qry
# read in the promoter information of the CMP genes

funcFindOverlapMotifBreakR=function(qry,sbj){
  k=findOverlaps(qry,sbj)
  match_hit=data.frame(names(qry)[queryHits(k)],names(sbj)[subjectHits(k)], stringsAsFactors=F)
  names(match_hit)= c('FeatureID','Variant')
  qrydf=as.data.frame(qry)
  sbjdf=as.data.frame(sbj)
  sbjdf=sbjdf[rownames(sbjdf) %in% match_hit$sbjID,]
  Strand=character(); GeneID=character();Transcription_Start=numeric()
  for (j in rownames(sbjdf)){
    ex=match_hit[match_hit$sbjID==j,"FeatureID"]
    gn=qrydf[rownames(qrydf)==ex,"GeneID"]
    ts=qrydf[rownames(qrydf)==ex,"TSS"]
    GeneID=append(GeneID,ex)
    Strand=append(Strand,st)
    Transcription_Start=append(Transcription_Start,ts)
  }
  sbjdf=data.frame(sbjdf,GeneID,Transcription_Start,Strand,stringsAsFactors=F)
  sbjdf=subset(sbjdf,select=c("GeneID","Transcription_Start","Strand"))
  return(sbjdf)
}

funcAnnotateMotifInfo=function(x,y){
  mat=matrix(NA,0,4); colnames(mat)=c("Location","Distance_from_TSS","Gene","TFs")
  if (is.element(x,y$varID)==F){
    out=c("",0,"","")
  } else {
    xy=y[y$varID==x,];z=xy[1,]
    TFs=unique(xy$geneSymbol);TFs=TFs[is.na(TFs)==F]
    TFs=paste(TFs,sep=", ",collapse = ", ")
    gen=z$NearbyGene
    ds=z$snpPos-z$NearbyTSS
    if(z$Gene_Strand=="+" & ds>=0){
      Loc="Downstream"
    } else if(z$Gene_Strand=="+" & ds<0){
      Loc="Upstream"
    } else if(z$Gene_Strand=="-" & ds<0){
      Loc="Downstream"
    } else if (z$Gene_Strand=="-" & ds>0){
      Loc="Upstream"
    } else {
      Loc="Exact"
    }
    out=c(Loc,ds,gen,TFs)
  }
  mat=rbind(mat,out)
  mat=data.frame(mat,row.names=NULL,stringsAsFactors=F)
  return(mat)
}

funcAnnotateMotifInfoWithoutLocation=function(x,y){
  mat=matrix(NA,0,3); colnames(mat)=c("Distance_from_TSS","Gene","TFs")
  if (is.element(x,y$varID)==F){
    out=c(0,"","")
  } else {
    xy=y[y$varID==x,];z=xy[1,]
    TFs=unique(xy$geneSymbol);TFs=TFs[is.na(TFs)==F]
    TFs=paste(TFs,sep=", ",collapse = ", ")
    gen=z$NearbyGene
    ds=z$snpPos-z$NearbyTSS
    out=c(ds,gen,TFs)
  }
  mat=rbind(mat,out)
  mat=data.frame(mat,row.names=NULL,stringsAsFactors=F)
  return(mat)
}

#GeneID Transcription_Start Strand                     UniqueID

# function to extract data for probands
# x is a list/characters of all proband IDs : x=c("TR-CMP-1","TR-CMP-2", etc)
# y is the dataframe containing info to extract

funcExtractProb=function(x,y){
  mat=matrix(NA,0,ncol(y))
  for (i in 1:nrow(y)){
    ydt=y[i,]
    pt=unlist(strsplit(ydt$SampleID,split=", "))
    prob=intersect(x,pt)
    if (length(prob)>=1){
      iprob=ydt
      iprob$noSamples=length(prob)
      iprob$SampleID=paste(prob,sep=", ",collapse = ", ")
      mat=rbind(mat,iprob)
    }
  }
  return(mat)
}


# Function to summarize NCV for probands
funcSummariseNCV=function(x,y){
  mat=matrix(NA,0,3); colnames(mat)=c("SampleID","nNCV","NCV")
  ydt=y[y$UniqueSample==x,]
  if(nrow(ydt)==0){
    a=0; b=""
    mat=rbind(mat,c(x,a,b))
  } else {
    a=nrow(ydt); b=paste(ydt$Gene,ydt$varID,sep="_"); b=paste(b,sep=", ",collapse = ", ")
    mat=rbind(mat,c(x,a,b))
  }
  mat=data.frame(mat,stringsAsFactors = F); rownames(mat)=NULL
  return(mat)
}

funcMotifEffect=function(x,y){
  refscore=y$scoreRef[x]; altscore=y$scoreAlt[x]
  if(refscore>altscore){
    ef="Motif_loss"
  } else {
    ef="Motif_gain"
  }
  return(ef)
}


funcCollapseTF=function(x){
  if(nrow(x)>=1){
    y=paste(unique(x$geneSymbol[is.na(x$geneSymbol)==F]),sep="|", collapse="|")
  } else {
    y=""
  }
  return(y)
}

# function that split TFs by strand and by effect on motifs
funcTFstrand=function(x){
  mat=matrix(NA,0,4)
  xplus=x[x$strand=="+",]
  xminus=x[x$strand=="-",]
  xgain=x[x$MotifEffect=="Motif_gain",]
  xloss=x[x$MotifEffect=="Motif_loss",]
  p=funcCollapseTF(xplus)
  m=funcCollapseTF(xminus)
  g=funcCollapseTF(xgain)
  l=funcCollapseTF(xloss)
  mat=rbind(mat,c(p,m,g,l)); colnames(mat)=c("PlusStrandTF","MinusStrandTF","motifGain","motifLoss")
  mat=data.frame(mat,stringsAsFactors = F); rownames(mat)=NULL
  return(mat)
}


funcReduceMBR=function(x,y){
  xdat=y[y$varID2==x,]
  if (nrow(xdat)>=1){
    adat=xdat[1,]
    adat$geneSymbol=paste(unique(xdat$geneSymbol[is.na(xdat$geneSymbol)==F]),sep="|",collapse = "|")
    adat$dataSource=paste(unique(xdat$dataSource),sep="|",collapse = "|")
  } else {
    adat=xdat
  }
  adat=adat[,c("geneSymbol","dataSource","varID2")];colnames(adat)[1]="TFs"
  st=funcTFstrand(xdat)
  adat=cbind(adat,st)
  return(adat)
}

# update CardiacTF
funcUpdateCardiacTF=function(x,y,z){
  tf=y[x,"TFs"]; tf=unlist(strsplit(tf,split="\\|"))
  zdt=read.delim(z,header=F,stringsAsFactors = F); colnames(zdt)="CardiacTF"
  tfCardiac=intersect(zdt$CardiacTF,tf)
  if(length(tfCardiac)>=1){
    tfCardiac=paste(tfCardiac,sep="|",collapse = "|")
  } else {
    tfCardiac=""
  }
  return(tfCardiac)
}


funcReturnValue=function(p,q){
  mat=matrix(NA,0,3)
  Score=c()
  q=q[p,]
  RegDB=q$RegulomeDB_Score
  if(length(grep(RegDB,c("1a","1b","1c","1d","1e","1f","2a","2b","2c")))>=1){
    a=1; Score=append(Score,"RegulomeDB")
  } else {
    a=0
  }
  motifBreakR=q$TFs
  if(motifBreakR=="" | is.na(motifBreakR)==T){
    b=0
  } else {
    b=1; Score=append(Score,"motifBreakR")
  }
  DeepSea=q$DeepSEA_pvalue
  if(DeepSea<=0.05){
    c=1; Score=append(Score, "DeepSea")
  } else {
    c=0
  }
  fathmm=q$FATHMM_MKL_NC
  if(fathmm>=0.5){
    d=1; Score=append(Score,"fathmm")
  } else {
    d=0
  }
  Score=paste(Score,sep=", ",collapse=", ")
  NC_Score=a+b+c+d;
  if(NC_Score>=3){
  Class="Tier1"
  } else {
  Class="Tier2"
  }
  mat=rbind(mat,c(NC_Score,Score,Class)); mat=data.frame(mat,stringsAsFactors = F); rownames(mat)=NULL
  colnames(mat)=c("NCV_score","NCV_tools","Class")
  return(mat)
}
