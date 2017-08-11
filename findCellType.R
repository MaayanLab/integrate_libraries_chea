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

## Getting cell lines (ie, disease descriptions) from online
for(i in 1:length(rownames(fill))) {
  cell_URL<- paste0('https://scicrunch.org/resources/Cell%20Lines/search?q=',fill[i,1],'&l=', fill[i,1], '&facet[]=Organism:Homo%20sapiens&facet[]=Organism:Mus%20musculus', collapse=NULL)
  b<-getURL(cell_URL,ssl.verifypeer=FALSE)
  write(b,"C:/Users/maayanlab1/Downloads/searchTest.txt")
  searchTestFile<-read_file("C:/Users/maayanlab1/Downloads/searchTest.txt")
  if(!grepl("We could not find",searchTestFile)) {
    disease<-strapplyc(b,"Disease:([^\"]*)")
    disease<-unlist(disease)
    fill[i,2]<-disease[1]
  }
}

