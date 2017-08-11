#a<-read.table("~/sig_files_7-6-17/MINT.sig",sep=" ",nrows=5473,stringsAsFactors = FALSE)

## Note: I edited DIP; KMT2A NA NA NA NA MLL/AF6 fusion is now KMT2A NA NA NA NA MLL/AF6_fusion
## Also edited IntAct; MLL/AF6 fusion NA NA NA NA KMT2A is now MLL/AF6_fusion NA NA NA NA KMT2A (NA NA NA NA 0 Binding 20153263)
## Also edited MINT; BCR/ABL fusion NA NA NA NA CTNNB1 is now BCR/ABL_fusion NA NA NA NA CTNNB1 (NA NA NA NA 0 Binding 17318191)


files<-list.files(path="~/sig_files_7-6-17") 

ppi<-data.frame(source_protein=character(),target_protein=character(),effect=character(),interaction=character(),PMID=character(),file=character(),stringsAsFactors = FALSE)

extractPPI<-function(extract,ppi,name) {
  #assign("ppi",ppi,.GlobalEnv)
  extract_temp<-data.frame(extract$V1,extract$V6,extract$V11,extract$V12,extract$V13,rep(name,length(rownames(extract))),stringsAsFactors = FALSE)
  #print(head(extract))
  colnames(extract_temp)=c("source_protein","target_protein","effect","interaction","PMID","file")
  return(extract_temp)
}

for(i in 1:length(files)) {
  print(basename(files[i]))
  if(grepl("2017",basename(files[i]))) {
    extract<-read.table(paste0("~/sig_files_7-6-17/",basename(files[i])),sep="\t",stringsAsFactors = FALSE)
    extract_temp<-extractPPI(extract,ppi,basename(files[i]))
    ppi<-rbind(ppi,extract_temp)
  }
  else {
    extract<-read.table(paste0("~/sig_files_7-6-17/",basename(files[i])),sep=" ",stringsAsFactors = FALSE)
    extract_temp<-extractPPI(extract,ppi,basename(files[i]))
    ppi<-rbind(ppi,extract_temp)
  }
  print(i)
}
## options(error=NULL)

## Extracting only tfs
ppiTfs<-ppi[(ppi[,1]%in%tfs)&(ppi[,2]%in%tfs),]

library(magrittr)

ppiTfs1 <- ppiTfs[!duplicated(ppiTfs[1:2]),]

letterSorted<-t(apply(ppiTfs2,1,sort))
letterSorted<-as.data.frame(ppiTfs3,stringsAsFactors = FALSE)

ppiTfs2<-ppiTfs1[!duplicated(letterSorted[1:2]),]

write.table(ppiTfs2,"~/ppiTfs.tsv",row.names = FALSE,col.names=TRUE)

