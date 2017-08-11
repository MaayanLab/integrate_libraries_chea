library(dplyr)
library(reshape2)

intersection<-unique(intersect(colnames(coexp),colnames(chea2))) 
intersection<-unique(intersect(intersection,colnames(encode2)))

p<-intersection%in%area4
p<-subset(intersection,p==T)

creedsIntersection<-creedsFisher[,gsub("-.*","",colnames(creedsFisher))%in%p==T]

## Subset datasets
datasetsListSmall<-list(cheaDash,encodeDash,coexp) #below for function changes colnames

datasetsListSmall1<-lapply(datasetsListSmall,function(dataset){ q<-(gsub("-.*","",colnames(as.data.frame(dataset)))%in%p)
                                                                print(q[1:5])                        
                                                                return(as.data.frame(dataset)[,q==T])})
#for combining datasets
for(i in 1:length(datasetsListSmall)){
  p<-gsub("-.*","",colnames(as.data.frame(datasetsListSmall1[[i]])))
  colnames(datasetsListSmall[[i]])<-p
}


d<-datasetsListSmall1[1]
d<-as.data.frame(d)

e<-datasetsListSmall3[2]
e<-as.data.frame(e)

f<-datasetsListSmall3[3]
f<-as.data.frame(f)

write.table(f,"~/datasetsListSmall3-Coexp.tsv",col.names=T)

## Find p-vals
pValSimple<-function(tf1,tf2) {
  
  tf1<-tf1[!is.na(tf1)]
  tf2<-tf2[!is.na(tf2)]
  
  contTable<-matrix(c((20000-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  return(pVal[1])
}

makepVals<-function(x) {
  pVals<-as.data.frame(do.call(cbind,(apply(creedsIntersection, 2, function(d) apply(x,2,pValSimple,tf1=d))))) 
  
  pVals<-apply(pVals,c(1,2),unlist)
  pVals<-as.data.frame(pVals,stringsAsFactors=F)
  return(pVals)
}

ptm<-proc.time()
datasetsListSmall2<-lapply(datasetsListSmall1,makepVals)
proc.time()-ptm

## Finding the minimum p-value of experiments with the same tfs
minimum<-function(a,b,dataset=dataset){ 
  rows<-grep(paste0("^",a,"$"),gsub("-.*","",rownames(as.data.frame(dataset))))
  print(rows)
  if(identical(rows,integer(0))) { # shouldn't happen right now
    return(1)
  }
  else {
    return(min(as.data.frame(dataset)[rows,b]))
  } 
}

consol<-function(dataset) {
  minPValsLarge<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset))), function(d) sapply(sort(unique(gsub("-.*","",rownames(as.data.frame(dataset))))),minimum,b=d,dataset=dataset))))) #
  colnames(minPValsLarge)<-colnames(dataset)
  return(minPValsLarge)
}

ptm<-proc.time()
datasetsListSmall3<-lapply(datasetsListSmall2,consol) # technically don't have to do this for coexp #FIX
proc.time()-ptm

## Rank for mean rank and top rank
datasetsListSmall4<-lapply(datasetsListSmall3,function(x) apply(x,2,rank))

u<-as.data.frame(datasetsList6[[1]])

## Do some method. Go to methods.R. Can use code below to save results.
datasetsList4<-list()

datasetsListMultipliedPValSmall<-combn(datasetsListSmall3,2,wrapperMethods,method=1)

datasetsListMeanRankSmall<-combn(datasetsListSmall4,2,wrapperMethods,method=2)
datasetsListTopRankSmall<-combn(datasetsListSmall4,2,wrapperMethods,method=3)

three<-method_df(as.data.frame(datasetsListMultipliedPValSmall[[1]]),as.data.frame(datasetsListSmall3[[3]]),method=1)
datasetsListMultipliedPValSmall[[length(datasetsListMultipliedPValSmall)+1]]<-three
datasetsList5<-list(as.data.frame(datasetsListSmall3[[1]]),as.data.frame(datasetsListSmall3[[2]]),as.data.frame(datasetsListSmall3[[3]]),as.data.frame(datasetsListMultipliedPValSmall[[1]]),as.data.frame(datasetsListMultipliedPValSmall[[2]]),as.data.frame(datasetsListMultipliedPValSmall[[3]]),as.data.frame(datasetsListMultipliedPValSmall[[4]]))

three<-minRankThree(as.data.frame(datasetsListSmall4[[1]]),as.data.frame(datasetsListSmall4[[2]]),as.data.frame(datasetsListSmall4[[3]]))
datasetsListMeanRankSmall[[length(datasetsListMeanRankSmall)+1]]<-three 
datasetsList5<-datasetsListMeanRankSmall
datasetsList5<-list(as.data.frame(datasetsListSmall4[[1]]),as.data.frame(datasetsListSmall4[[2]]),as.data.frame(datasetsListSmall4[[3]]),as.data.frame(datasetsListMeanRankSmall[[1]]),as.data.frame(datasetsListMeanRankSmall[[2]]),as.data.frame(datasetsListMeanRankSmall[[3]]),as.data.frame(datasetsListMeanRankSmall[[4]]))

three<-method_df(as.data.frame(datasetsListTopRankSmall[[1]]),as.data.frame(datasetsListSmall4[[3]]),method=3)
datasetsListTopRankSmall[[length(datasetsListTopRankSmall)+1]]<-three
datasetsList5<-datasetsListTopRank
datasetsList5<-list(as.data.frame(datasetsListSmall4[[1]]),as.data.frame(datasetsListSmall4[[2]]),as.data.frame(datasetsListSmall4[[3]]),as.data.frame(datasetsListTopRankSmall[[1]]),as.data.frame(datasetsListTopRankSmall[[2]]),as.data.frame(datasetsListTopRankSmall[[3]]),as.data.frame(datasetsListTopRankSmall[[4]]))

m<-as.data.frame(datasetsListSmallCombine[[4]])
n<-as.data.frame(datasetsListMultipliedPVal[[1]])
o<-as.data.frame(datasetsListTopRank[[4]])

datasetsListMultipliedPValSmall<-combn(datasetsListSmall3,2,wrapperMethods,method=1)

#
temp<-list()
datasetsListCombine2<-list()
sample<-combn(datasetsListSmall,2,wrapperIntersect)
datasetsListSmallCombine<-sample

threeMax<-findMaxIntersect(as.data.frame(datasetsListSmallCombine[[1]]),as.data.frame(datasetsListSmall[[3]]))
threeCombine<-make_intersect_df(as.data.frame(datasetsListSmallCombine[[1]]),as.data.frame(datasetsListSmall[[3]]),threeMax)
datasetsListSmallCombine[[length(datasetsListSmallCombine)+1]]<-threeCombine #

datasetsListSmallCombine2<-lapply(datasetsListSmallCombine,makepVals)
datasetsListSmallCombine3<-list(as.data.frame(datasetsListSmall3[[1]]),as.data.frame(datasetsListSmall3[[2]]),as.data.frame(datasetsListSmall3[[3]]),as.data.frame(datasetsListSmallCombine2[[1]]),as.data.frame(datasetsListSmallCombine2[[2]]),as.data.frame(datasetsListSmallCombine2[[3]]),as.data.frame(datasetsListSmallCombine2[[4]]))
datasetsList5<-datasetsListSmallCombine3

## Combine
ranks<-as.data.frame(do.call(cbind,(lapply(datasetsList6, function(t) t(sapply(1:length(colnames(as.data.frame(t))),ranked,d=t)))))) #length(colnames(as.data.frame(t)))
rownames(ranks)=colnames(creedsIntersection)
colnames(ranks)=c("ChEA","ChEA Control","ENCODE","ENCODE Control","Coexp","Coexp Control","ChEA Encode","ChEA ENCODE Control","ChEA Coexp","ChEA Coexp Control","ENCODE Coexp","ENCODE Coexp Control","All three","All three Control")
ranks<-apply(ranks,2,as.numeric) 
ranks<-as.data.frame(ranks,stringsAsFactors = F)

write.table(ranks,"~/ranksSmallIntersectionScaled.tsv",col.names=T)

ranked<-function(x,d){ 
  a<-which(rownames(as.data.frame(d))==gsub("-.*","",colnames(as.data.frame(d))[x]))
  if(identical(a,integer(0))) {
    r<-NA
  }
  else {
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  print(paste(r,s))
  return(data.frame(exp=r,control=s))
}

## plot
library(reshape2)
ranksMeltNormalize<-melt(ranks)

ggplot(data=subset(ranksMeltNormalize,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,3.5),expand = c(0,0))
