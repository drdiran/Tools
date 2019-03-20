#!/usr/bin/env Rscript

library(getopt)
library(XML)
library(httr)
library(xml2)

### FUNCTIONS ######################################################################################
usage <- function(){
  usage.value <-
    '
  RegulomeDB_Query_Pipeline.R
  -help/-h (help)
  -RegDB_File/-r (Output NCV_pipeline)
  ';
  cat(usage.value);
}

params <- matrix(
  c(
    'help', 'h', 0, 'logical',
    'RegDB_File', 'r', 1, 'character'
  ),
  byrow = TRUE,
  ncol = 4
);

opt <- getopt(params);

if(is.null(opt$RegDB_File)){
  print("no input file provided");
  usage();
  stop(options("show.error.messages" = FALSE));
}


# regulomedb url for web scrapping
url="http://www.regulomedb.org/snp"

snpCord=opt$RegDB_File
fname=unlist(strsplit(snpCord,split="_Input"))[1]
OutputName=paste(fname,"RegDB_Scores.RData",sep="_")
OutputName_2=paste(fname,"RegDB_Scores.txt",sep="_")
#webqry="RegulomeDB_Input_For_Web_Qry.bed"


# 
# funcQueryRegulomeDB=function(url,x){
#   x=gsub(":","/",x)
#   qryURL=paste(url,x,sep="/")
#   download=try(htmlParse(qryURL))
#   if(class(download)[1]=="try-error"){
#     score="Invalid_coordinate"
#   } else {
#     score=xpathSApply(download, "//h2", xmlValue)[1]
#     score=gsub("Score: \n","",score)
#   }
#   return(score)
# }


funcQueryRegulomeDB_2=function(url,x){
  x=gsub(":","/",x)
  qryURL=paste(url,x,sep="/")
  download=tryCatch( htmlParse(qryURL), XMLError = function(e){msg=e$message})
  if(class(download)[1]=="character"){
    score="Invalid_coordinate"
  } else {
    score=xpathSApply(download, "//h2", xmlValue)[1]
    score=gsub("Score: \n","",score)
  }
  return(score)
}

inPut=read.delim(snpCord,header=F,sep="\t",stringsAsFactors = F)
inPut$RegDBinput=paste(inPut$V1,inPut$V2,sep="/")

if(length(inPut$RegDBinput[grepl("^chr",inPut$RegDBinput)==T])==0){
  inPut$RegDBinput=paste("chr",inPut$RegDBinput,sep="")
}

inPut=inPut[!duplicated(inPut$RegDBinput),]

inPut$RegulomeDB_Score=do.call(rbind,lapply(inPut$RegDBinput,funcQueryRegulomeDB_2,url=url))
RegScore=inPut
write.table(RegScore,OutputName_2,row.names=F,sep="\t",quote=F)
save(RegScore,file=OutputName)


