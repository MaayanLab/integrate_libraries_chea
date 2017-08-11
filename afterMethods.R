## rank again
datasetsList5<-lapply(datasetsList5,function(x) apply(x,2,rank))

j<-as.data.frame(datasetsList5[[5]])
k<-as.data.frame(datasetsList6[[5]])
l<-as.data.frame(datasetsList5[[7]])

## normalize/scale

datasetsList6<-lapply(datasetsList5,function(a)apply(a,2,oldNormalize))

oldNormalize<-function(x) {
  largest<-max(x)
  return(sapply(x,function(y) return(y/largest)))
}

# rank
ranks<-as.data.frame(do.call(cbind,(lapply(datasetsList6, function(t) t(sapply(1:length(colnames(as.data.frame(t))),ranked,d=t)))))) #length(colnames(as.data.frame(t)))
rownames(ranks)=colnames(creedsFisher)
colnames(ranks)=c("ChEA","ChEA Control","ENCODE","ENCODE Control","Coexp","Coexp Control","ChEA Encode","ChEA ENCODE Control","ChEA Coexp","ChEA Coexp Control","ENCODE Coexp","ENCODE Coexp Control","All three","All three Control")
ranks<-apply(ranks,2,as.numeric)
ranks<-as.data.frame(ranks,stringsAsFactors = F)
ranks<-as.matrix(ranks)

write.table(ranks,"~/ranksTopRankScaled.tsv",col.names = T,row.names=T)
ranksSmallIntersectionScaled<-ranks

ranked<-function(x,d){ 
  a<-which(rownames(as.data.frame(d))==gsub("-.*","",colnames(as.data.frame(d))[x]))
  if(identical(a,integer(0))) {
    r<-NA
  }
  else {
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  return(data.frame(exp=r,control=s))
}

#change order. Doesn't actually change the graph for some reason.
# allThreeControl<-ranksMeanRankPlotScaled %>% filter(variable=="All three Control")
# ranksMeanRankPlotScaled1<-ranksMeanRankPlotScaled %>% filter(variable!="All three Control")
# ranksMeanRankPlotScaled1<-rbind(ranksMeanRankPlotScaled1,allThreeControl)

#plot
# ranks=ranksTopRankScaled[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]
# ranks<-ranksMeanRankScaled
ranks<-ranksSmallIntersectionScaled
ranksMeltNormalize<-melt(ranks)
write.table(ranksMeltNormalize,"~/ranksSmallIntersectionPlotScaled.tsv",sep="\t",col.names = T)
ranksSmallIntersectionPlotScaled<-ranksMeltNormalize

ggplot(data=subset(ranksMeanRankPlotScaled,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,3.5),expand = c(0,0))#+scale_colour_manual(values=c("white","white","white","white","white","blue","white","white","white","white","white","white","white","white")) 

ranks_subset<- ranksSmallIntersectionPlotScaled %>% filter(variable %in% c("ENCODE Coexp","ChEA Coexp")==T)

barplot(as.matrix(ranks_subset),main="Ranks",xlab="Rank of TF",legend=colnames(ranks))
barplot(as.matrix(ranks_subset),main="Ranks",beside=T)
## 

test<-as.data.frame(do.call(cbind,(apply(ranksTopRankUnscaled,2,function(x) sapply(1:2,function(t) {b<-(sum(unlist(x)==t,na.rm=T))
print(unlist(x)[1:5])
return(b)})))))

example <-function(x) {
  
  print(unlist(x)[1:5])
  sapply(1:417,function(a) return(sum(unlist(x)==a,na.rm=T)))
}
