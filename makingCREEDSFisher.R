install.packages("stats")

library(stats)
library(stringi)
library(dplyr)

## Find out how many rows the data frame should be
b<-apply(creeds,1,greatest)
max(b)

greatest<-function(x){
  up<-stri_split_fixed(x[3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(x[4],",")
  down<-sapply(down,toupper)
  
  return(length(up)+length(down))
}
## Format creeds to something that looks like chea and encode

format<-function(d){
  up<-stri_split_fixed(d[3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(d[4],",")
  down<-sapply(down,toupper)
  name<-paste(d[2],d[1],d[5],d[6],collapse=":")
  temp<-data.frame(name=character(601),stringsAsFactors=FALSE)
  temp[,1]<-NA
  temp[1:length(up),1]<-up
  temp[(length(up)+1):(length(up)+length(down)),1]<-down
  colnames(temp)=name
  #print(typeof(temp))
  return(temp)
}

creedsFisher<-as.data.frame(do.call(cbind,(apply(creeds,1, function(d) format(d)))))
colnames(creedsFisher)<-gsub(" ","-",colnames(creedsFisher))

colnames(creedsFisher)<-colnames(as.data.frame(datasetsList4[[1]]))

duplicates<-colnames(creedsFisher)[duplicated(colnames(creedsFisher))==T]

for(i in 1:length(duplicates)) { #:length(duplicates)
  a<-grep(duplicates[i],colnames(creedsFisher))
  for(i in 1:length(a)) {
    colnames(creedsFisher)[a[i]]<-paste(colnames(creedsFisher)[a[i]],i,sep = "-") #colnames(creedsFisher)[a[i]]
  }
}

## Lazy way of accounting for "VDR-GSE2421-Bone-marrow-derived-macrophages-(BMDMs)-KO" because parentheses
colnames(creedsFisher)[119]="VDR-GSE2421-Bone-marrow-derived-macrophages-(BMDMs)-KO-1"
colnames(creedsFisher)[120]="VDR-GSE2421-Bone-marrow-derived-macrophages-(BMDMs)-KO-2"

write.table(creedsFisher,"~/creedsFisher.tsv",col.names=T,sep="\t")