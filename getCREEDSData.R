# Need: tfs, allCellTypes, dataset_ids

# install.packages("httr")
# install.packages("jsonlite")
# install.packages("stringi")

## Get IDs
library(httr)
library(jsonlite)
library(stringi)

CREEDS_URL <- 'http://amp.pharm.mssm.edu/CREEDS/'

ids <- character()
count=0
tf<-0

for(i in 1:length(tfs)) {
  tf<-tf+1
  print(tf)
  
  r <- GET(paste0(CREEDS_URL, 'search'), query=list(q = tfs[i]))
  response <- fromJSON(httr::content(r, 'text'))
  if(length(response)!=0 && "hs_gene_symbol" %in% colnames(response)) {
    response = subset(response,response$hs_gene_symbol==tfs[i])
    response<-subset(response,!(is.na(response["hs_gene_symbol"])))
    if(length(response$id)!=0) {
      ids<-append(ids, response$id)
      count<-count+1
    }
  }
}

## Make table
creeds<-matrix(NA,count,6)
creeds<-as.data.frame(creeds)

count1=0
## Put into table
for(i in 1:length(ids)) {
  count1<-count1+1
  print(count1)
  r <- GET(paste0(CREEDS_URL, 'api'), query=list(id= ids[i]))
  sig <- fromJSON(httr::content(r, 'text'))
  
  # numeric<-as.numeric(sig$up_genes)
  # indices<-which(is.na(numeric))
  # 
  # numeric1<-as.numeric(sig$down_genes)
  # indices1<-which(is.na(numeric1))
  
  creeds[i,1]<-sig$geo_id
  creeds[i,2]<-sig$hs_gene_symbol
  creeds[i,3]<-paste0(sig$up_genes[indices],collapse=",")
  creeds[i,4]<-paste0(sig$down_genes[indices1],collapse=",")
  # creeds[i,3]<-paste0(sig$up_genes[,1],collapse=",") doesn't work because sometimes there isn't any up genes, so no dimensions.
  # creeds[i,4]<-paste0(sig$down_genes[,1],collapse=",")
  creeds[i,5]<-sig$cell_type
  creeds[i,6]<-sig$pert_type
  
}

## Making creedsSep, with separate interactions for a tf and all targets
creedsSep<-data.frame(a=character(), b=character(), c=character(), d=character(),e=character(), f=character(), g=character(), h=character(), i=character(), stringsAsFactors = FALSE)
colnames(creedsSep) = c("source","interaction","target","experiment (GEO series or PMID)","cell line","cell type","perturbation","organism","source:target")


makeCreedsSep<-function(d) { #:length(rownames(creeds))
  print(d[2])
  up<-stri_split_fixed(d[3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(d[4],",")
  down<-sapply(down,toupper)
  
  sourcesUp = rep(d[2],length(up))
  interactionsUp = rep("up",length(up))
  tempUp<-data.frame(sourcesUp,interactionsUp,up)
  colnames(tempUp) = c("a","b","c")
  
  sourcesDown = rep(d[2],length(down))
  interactionsDown = rep("down",length(down))
  tempDown<-data.frame(sourcesDown,interactionsDown,down)
  colnames(tempDown) = c("a","b","c")
  
  creedsSepTemp = rbind(tempUp,tempDown)
  
  creedsSepTemp[,4:9]<-NA
  creedsSepTemp[,4]<-d[1]
  creedsSepTemp[,5]<-d[5]
  #creedsSepTemp[,6]<-cellTypes2[which(creedsSepTemp[i,5]==cellTypes2[,1]),2]
  creedsSepTemp[,7]<-d[6]
  #creedsSepTemp[i,9]<-paste0(creedsSepTemp[i,1],":",creedsSepTemp[i,3],collapse=NULL) # may see problems here # forgot [i,9]-should be fine now.
  colnames(creedsSepTemp) = c("a","b","c","d","e","f","g","h","i")
  return(creedsSepTemp)
  
}
creedsSep<-as.data.frame(do.call(rbind,(apply(creeds,1, function(d) makeCreedsSep(d)))))
colnames(creedsSep)=c("source","interaction","target","experiment (GEO series or PMID)","cell line","cell type","perturbation","organism","source:target")

## Getting cell types (added one, difference between creedsSep and its children)
allCellTypes<-read.table("C:/Users/maayanlab1/Downloads/allCellTypes - allCellTypes.tsv",sep="\t",stringsAsFactors = FALSE) # From Google Sheets
allCellTypes[719,1]<-"strain/background: C57BL/6 genotype/variation: Zbtb7a Flox/Flox Mx1-Cre+ (LRF homozygous knockout) age: 4 weeks cell type: FACS-sorted LT-HSCs (LSK IL7Ra-Flt3-CD150+CD48-) time point: 9 days after 1st pIpC injection Treatment protocol\tLT-HSCs were FACS-sorted from three Lrf knockout (Zbtb7a Flox/Flox Mx1-Cre+) and two control (Zbtb7a Flox/+ Mx1-Cre+) mice, nine days after the first pIpC injection. pIpC (250&Icirc;&frac14;g) was injected twice to induce Cre recombinase."
# \t was missing when it imported.

a<-sapply(creedsSep[,5],function(x){return(allCellTypes[which(x==allCellTypes[,1]),2])} )
creedsSep[,6]<-a

## Getting whether it's mouse or human 
dataset_ids<-read.table("C:/Users/maayanlab1/Downloads/dataset_ids - dataset_ids.tsv",sep="\t",stringsAsFactors = FALSE) # From Google Sheets
dataset_ids<-dataset_ids[!duplicated(dataset_ids[,1]),]
creedsSep[creedsSep[,4]=="GSE2715",8]<-"Homo sapiens" #There is GSE27157,GSE2715,GSE27159 and grep identifies all of them.

ptm <- proc.time()
for(i in 1:length(rownames(creedsSep))) { 
  #print(i)
  num<-grep(creedsSep[i,4],dataset_ids[,1],fixed=TRUE)
  if(length(num)) {
    creedsSep[i,8]<-dataset_ids[num[1],2]
  }
}
proc.time() - ptm

## Fix source:target (small problem in code above)
a<-apply(creedsSep,1,function(x){return(paste0(x[1],":",x[3],collapse=NULL))} )
creedsSep[,9]<-a

write.table(creedsSep,"~/creedsSepFull.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
creedsSep<-read.table("~/creedsSepFull.tsv",header=TRUE,sep="\t",stringsAsFactors = FALSE)

##creeds, just tf-tf interactions #should fix if needed
library(dplyr)
creedsTfs<-creedsSep %>% 
           filter(V3%in%tfs) %>%
           rename(.,source=V1,interaction=V2,target=V3,"experiment (GEO series or PMID)"=V4,"cell line"=V5, "cell type"=V6, "perturbation"=V7, "organism"=V8, "source:target"=V9)
write.table(.,"~/creedsTfs.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
