install.packages("reshape2")
install.packages("rlist")

library(dplyr)
library(rlist)

## Get p-values for the creeds experiment with the tfs (in the intersection) from the dataset

#genes<-unique(union(unlist(datasetIntersection),unlist(creedsFisher)))%>%na.omit() set as 20,000

datasetsList<- list(cheaDash,encodeDash,coexp)

makepVals<-function(x,which_creeds) {
  
  pValSimple<-function(tf1,tf2) {
    
    tf1<-tf1[!is.na(tf1)]
    tf2<-tf2[!is.na(tf2)]
    # print(paste(head(tf1),head(tf2),sep="---"))
    
    contTable<-matrix(c((20000-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
    pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
    return(pVal[1])
  }
  
  pVals<-as.data.frame(do.call(cbind,(apply(which_creeds, 2, function(d) apply(x,2,pValSimple,tf1=d))))) 
  
  pVals<-apply(pVals,c(1,2),unlist)
  pVals<-as.data.frame(pVals,stringsAsFactors=F)
  return(pVals)
}

ptm<-proc.time()
datasetsList1<-lapply(datasetsList,makepVals,creedsFisher)
proc.time()-ptm

## trying out
d<-datasetsList1[1]
d<-as.data.frame(d)

e<-datasetsList1[2]
e<-as.data.frame(e)

f<-datasetsList1[3]
f<-as.data.frame(f)

## Finding the minimum p-value of experiments with the same tfs

consol<-function(dataset) {
  
  minimum<-function(a,b,dataset=dataset){ 
    rows<-grep(paste0("^",a,"$"),gsub("-.*","",rownames(as.data.frame(dataset))),ignore.case = T)
    
    if(identical(rows,integer(0))) { # shouldn't happen right now
      return(1)
    }
    else {
      return(min(as.data.frame(dataset)[rows,b]))
    } 
  }
  
  minPValsLarge<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset))), function(d) sapply(sort(unique(gsub("-.*","",rownames(as.data.frame(dataset))))),minimum,b=d,dataset=dataset))))) #
  colnames(minPValsLarge)<-colnames(dataset)
  return(minPValsLarge)
}

ptm<-proc.time()
datasetsList2<-lapply(datasetsList1,consol) # technically don't have to do this for coexp
proc.time()-ptm

## Try out
a<-datasetsList2[1]
a<-as.data.frame(a)

b<-datasetsList2[2]
b<-as.data.frame(b)

c<-datasetsList2[3]
c<-as.data.frame(c)

# write.table(c,"~/datasetsList2-Coexp.tsv",col.names=T)

## To change column names
for(i in 1:length(datasetsList2)){
  colnames(datasetsList2[[i]])<-colnames(creedsFisher)
}

## Rank
datasetsList3<-lapply(datasetsList2,function(x) apply(x,2,rank))

g<-datasetsList3[[1]]
g<-as.data.frame(g)

h<-datasetsList3[2]
h<-as.data.frame(h)

i<-datasetsList3[3]
i<-as.data.frame(i)

## Combine datasets- see methods.R. Set result to datasetsList5

## rank again
datasetsList5<-lapply(datasetsList5,function(x) apply(x,2,rank))

j<-as.data.frame(datasetsList0_ppi_new_creeds[[1]])
k<-as.data.frame(datasetsList6[[5]])
l<-as.data.frame(datasetsList5[[7]])

## normalize/scale

datasetsList6<-lapply(datasetsList5,function(a)apply(a,2,oldNormalize))

oldNormalize<-function(x) {
  largest<-max(x)
  return(sapply(x,function(y) return(y/largest)))
}

## rank
ranks<-as.data.frame(do.call(cbind,(lapply(datasetsList6, function(t) t(sapply(1:length(colnames(as.data.frame(t))),ranked,d=t)))))) #length(colnames(as.data.frame(t)))
rownames(ranks)=colnames(creedsFisher)
colnames(ranks)=c("ChEA","ChEA Control","ENCODE","ENCODE Control","Coexp","Coexp Control","ChEA Encode","ChEA ENCODE Control","ChEA Coexp","ChEA Coexp Control","ENCODE Coexp","ENCODE Coexp Control","All three","All three Control")
ranks<-apply(ranks,2,as.numeric)
ranks<-as.data.frame(ranks,stringsAsFactors = F)
ranks<-as.matrix(ranks)

write.table(ranks,"~/ranksTopRankScaled.tsv",col.names = T,row.names=T)
ranksSmallIntersectionScaled<-ranks

ranked<-function(x,d){ 
  b<-gsub("-.*","",colnames(as.data.frame(d))[x])
  a<-grep(paste0("^",b,"$"),rownames(as.data.frame(d)),ignore.case = T) #check datasetsList0
  #a<-which(rownames(as.data.frame(d))==gsub("-.*","",colnames(as.data.frame(d))[x]))
  if(identical(a,integer(0))) {
    r<-NA
  }
  else {
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  return(data.frame(exp=r,control=s))
}

## plot
maximum<-max(as.vector(unlist(ranks)),na.rm = T)
ranksNormalize<-apply(ranks,c(1,2),normalize,maximum)
ranksNormalize<-as.data.frame(ranksNormalize,stringsAsFactors = F)

normalize<-function(x,maximum) {
  x<-as.numeric(x)
  return(x/maximum)
}

ranksMeltNormalize<-melt(ranksNormalize)

ggplot(data=subset(ranksMeltNormalize,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,3.5),expand = c(0,0))

