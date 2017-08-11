# Need: technically noOrganism, but you can just use all datasets needed (may take long time though)

#install.packages("XML")
install.packages("RCurl")
#install.packages("xml2")
install.packages("readr")

library(httr)
library(jsonlite)
#library(XML)
library(RCurl)
#library(xml2)
library(readr)
library(gsubfn)

dataset_ids<-unique(noOrganism[,4])
dataset_ids<-as.data.frame(dataset_ids,stringsAsFactors = FALSE)
dataset_ids[,2]<-character(0)

for(i in 1:length(rownames(dataset_ids))) { ## this took 12.5 minutes. use the file.
  GEO_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=('
  b<-getURL(paste0(GEO_URL,dataset_ids[i,1],')&retmax=10&usehistory=y', collapse=NULL),ssl.verifypeer=FALSE)
  query_key<-strapplyc(b,"<QueryKey>(.*)</QueryKey>")
  web_env<-strapplyc(b,"<WebEnv>(.*)</WebEnv>")
  
  second_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key='
  
  c<-getURL(paste0(second_URL, query_key,'&WebEnv=',web_env, collapse=NULL),ssl.verifypeer=FALSE)
  organism<-strapplyc(c,"Organism:\t(.*?)\nType:")
  organism<-unlist(organism)
  organism<-organism[1]
  dataset_ids[i,2]<-organism
  print(i)
}

write.table(dataset_ids,"C:/Users/maayanlab1/Downloads/dataset_ids.tsv",sep="\t",row.names = FALSE,col.names=FALSE)

dataset_ids<-read.table("C:/Users/maayanlab1/Downloads/dataset_ids.tsv",sep="\t",stringsAsFactors = FALSE)
#colnames(dataset_ids)=c("a","b")

dataset_ids[dataset_ids[,2]=="Homo sapiens",2]<-"human"
dataset_ids[dataset_ids[,2]=="Mus musculus",2]<-"mouse"

