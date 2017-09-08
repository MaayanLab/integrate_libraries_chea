## Load files

downloaded_file<-read.ttsv("~/Single_Gene_Perturbations_from_GEO_up.txt")
downloaded_file<-as.data.frame(downloaded_file)
downloaded_file = downloaded_file[2:length(rownames(downloaded_file)),]
downloaded_file = as.data.frame(gsub(",1.0","",as.matrix(downloaded_file)),stringsAsFactors = FALSE)

new_creeds_up<-downloaded_file

write.table(downloaded_file,"~/bioplex.tsv",sep="\t",col.names = T)

## Install SuperExactTest
install.packages("SuperExactTest")
library(SuperExactTest)

## Make new_creeds
library(dplyr)
new_creeds<-full_join(new_creeds_up,new_creeds_down)
a<-gsub("\\.","-",colnames(new_creeds))
a<-toupper(colnames(new_creeds))
colnames(new_creeds)<-a

install.packages("reshape2")
install.packages("rlist")

library(dplyr)
library(rlist)

## Get p-values 

datasetsList_ppi<- list(bioplex,humap) # The p-values for the other datasets are in datasetsList1

datasetsList0_ppi_new_creeds<-lapply(datasetsList,makepVals,new_creeds)
datasetsList1_ppi_new_creeds<-list(datasetsList0_ppi_new_creeds[[1]],datasetsList0_ppi_new_creeds[[2]],datasetsList0_ppi_new_creeds[[3]],bioplex_pval,humap_pval)

## trying out
d<-datasetsList0_ppi_new_creeds[[3]]
d<-as.data.frame(d)

## Adding in the other datasets
datasetsList2_ppi<-list(datasetsList2[[1]],datasetsList2[[2]],datasetsList2[[3]],datasetsList1_ppi[[1]],datasetsList1_ppi[[2]])

## Finding the minimum p-value of experiments with the same tfs. Even though HuMap and Bioplex don't have experiment names, 
## should still do. This could be changed; don't necessarily have to take the minimum p-value.

ptm<-proc.time()
datasetsList2_ppi<-lapply(datasetsList1_ppi,consol) 
proc.time()-ptm

ptm<-proc.time()
datasetsList2_ppi_new_creeds<-lapply(datasetsList1_ppi_new_creeds,consol) 
proc.time()-ptm

# write.table(c,"~/datasetsList2-Coexp.tsv",col.names=T)

## To change column names to uppercase
for(i in 1:length(datasetsList2_ppi_new_creeds)){
  colnames(datasetsList2_ppi_new_creeds[[i]])<-toupper(colnames(datasetsList2_ppi_new_creeds[[i]]))
}

## Rank
datasetsList3_ppi<-lapply(datasetsList2_ppi,function(x) apply(x,2,rank))
datasetsList3_ppi_new_creeds<-lapply(datasetsList2_ppi_new_creeds,function(x) apply(x,2,rank))

g<-datasetsList3_ppi[[1]]
g<-as.data.frame(g)

h<-datasetsList3_ppi[2]
h<-as.data.frame(h)

i<-datasetsList3_ppi[3]
i<-as.data.frame(i)

## Combine datasets- see revised_methods.R. Set result to datasetsList5

datasetsList5_ppi<-datasetsListMultipliedPVal_ppi

## Rank again
datasetsList5_ppi<-lapply(datasetsListMultipliedPVal_ppi_2_ranks,function(x) apply(x,2,rank))
datasetsList5_ppi_new_creeds<-lapply(datasetsList5_ppi_new_creeds,function(x) apply(x,2,rank))

j<-as.data.frame(datasetsListMultipliedRanks_ppi_2[[1]])
k<-as.data.frame(datasetsList2_ppi_new_creeds[[1]])
l<-as.data.frame(datasetsList6_ppi_new_creeds[[3]])

## normalize/scale

datasetsList6_ppi<-lapply(datasetsList3_ppi,function(a)apply(a,2,oldNormalize))
datasetsList6_ppi_new_creeds<-lapply(datasetsList3_ppi_new_creeds,function(a)apply(a,2,oldNormalize))

## Assemble ranks for all datasets into one data frame
ptm<-proc.time()
ranks<-as.data.frame(do.call(cbind,(lapply(datasetsList6_ppi_new_creeds, function(t) t(sapply(1:length(colnames(as.data.frame(t))),ranked,d=t)))))) #length(colnames(as.data.frame(t)))
proc.time()-ptm

ranked<-function(x,d){ 
  b<-gsub("-.*","",colnames(as.data.frame(d))[x])
  a<-grep(paste0("^",b,"$"),rownames(as.data.frame(d)),ignore.case = T)
  if(identical(a,integer(0))) {
    r<-NA
  }
  else {
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  return(data.frame(exp=r,control=s))
}

rownames(ranks)=colnames(new_creeds)
colnames(ranks)=c("ChEA","ChEA Control","ENCODE","ENCODE Control","Coexp","Coexp Control","Bioplex","Bioplex Control","HuMap","HuMap Control")
ranks<-apply(ranks,2,as.numeric)
ranks<-as.data.frame(ranks,stringsAsFactors = F)
ranks<-as.matrix(ranks)

ranks_ppi_new_creeds<-ranks
write.table(ranks,"~/ranks_ppi_new_creeds.tsv",col.names = T,row.names=T)

## plot
maximum<-max(as.vector(unlist(ranks)),na.rm = T)
ranksNormalize<-apply(ranks,c(1,2),normalize,maximum)
ranksNormalize<-as.data.frame(ranksNormalize,stringsAsFactors = F)

normalize<-function(x,maximum) {
  x<-as.numeric(x)
  return(x/maximum)
}

ranksMeltNormalize<-melt(ranksNormalize)

ggplot(data=subset(ranksMeltNormalize,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + labs(color="Datasets") + labs(x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,3.8),expand = c(0,0)) +
  scale_colour_manual(values = c("deepskyblue","deepskyblue","seagreen","seagreen","yellow","yellow","navy","navy","red3","red3"))

