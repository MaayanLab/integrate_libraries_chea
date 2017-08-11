## 1. multiplied p-value 2. mean rank 3. top rank 

wrapperMethods<-function(list,method) {
  dataset1<-as.data.frame(list[[1]])
  dataset2<-as.data.frame(list[[2]])
  toAdd<-method_df(dataset1,dataset2,method)
  datasetsList4[[length(datasetsList4)+1]]<-toAdd
  return(datasetsList4)
}

method_df<-function(dataset1,dataset2,method) {
  
  rows1<-gsub("-.*","",rownames(as.data.frame(dataset1)))
  rows2<-gsub("-.*","",rownames(as.data.frame(dataset2)))
  together<-as.vector(sort(union(rows1,rows2)))
  
  test<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset1))),function(d) sapply(1:length(together),method_compute,b=d,dataset1=dataset1,dataset2=dataset2,together=together,method=method))))) #length(colnames(as.data.frame(dataset1))) #length(rownames(as.data.frame(newDataset)))
  colnames(test)<-colnames(dataset1)
  rownames(test)<-together
  return(test)
}

method_compute <-function(a,b,dataset1,dataset2,together,method) {
  c<-paste0("^",together[a],"$")
  first<-grep(c,rownames(as.data.frame(dataset1)))
  second<-grep(c,rownames(as.data.frame(dataset2)))
  
  if(identical(first,integer(0))&&!identical(second,integer(0))) {
    return((as.data.frame(dataset2))[second,b])
  } else if(identical(second,integer(0))&&!identical(first,integer(0))) {
    return((as.data.frame(dataset1))[first,b])
  } else {
    if(identical(method,1)){
      return((as.data.frame(dataset1)[first,b])*(as.data.frame(dataset2)[second,b]))
    }
    if(identical(method,2)){
      return((as.data.frame(dataset1)[first,b]+as.data.frame(dataset2)[second,b])/2)
    }
    if(identical(method,3)){
      return(min(as.data.frame(dataset1)[first,b],as.data.frame(dataset2)[second,b]))
    }
  }
  
}

datasetsList4<-list()

ptm<-proc.time()
datasetsListMultipliedPVal<-combn(datasetsList2,2,wrapperMethods,method=1)
proc.time()-ptm

datasetsListMeanRank<-combn(datasetsList3,2,wrapperMethods,method=2)
datasetsListTopRank<-combn(datasetsList3,2,wrapperMethods,method=3)

three<-method_df(as.data.frame(datasetsListMultipliedPVal[[1]]),as.data.frame(datasetsList2[[3]]),method=1)
datasetsListMultipliedPVal[[length(datasetsListMultipliedPVal)+1]]<-three
datasetsList5<-list(as.data.frame(datasetsList2[[1]]),as.data.frame(datasetsList2[[2]]),as.data.frame(datasetsList2[[3]]),as.data.frame(datasetsListMultipliedPVal[[1]]),as.data.frame(datasetsListMultipliedPVal[[2]]),as.data.frame(datasetsListMultipliedPVal[[3]]),as.data.frame(datasetsListMultipliedPVal[[4]]))

three<-meanRankThreeWrapper(as.data.frame(datasetsList3[[1]]),as.data.frame(datasetsList3[[2]]),as.data.frame(datasetsList3[[3]]))
datasetsListMeanRank[[length(datasetsListMeanRank)+1]]<-three 
datasetsList5<-datasetsListMeanRank
datasetsList5<-list(as.data.frame(datasetsList3[[1]]),as.data.frame(datasetsList3[[2]]),as.data.frame(datasetsList3[[3]]),as.data.frame(datasetsListMeanRank[[1]]),as.data.frame(datasetsListMeanRank[[2]]),as.data.frame(datasetsListMeanRank[[3]]),as.data.frame(datasetsListMeanRank[[4]]))

three<-method_df(as.data.frame(datasetsListTopRank[[1]]),as.data.frame(datasetsList3[[3]]),method=3)
datasetsListTopRank[[length(datasetsListTopRank)+1]]<-three
datasetsList5<-datasetsListTopRank
datasetsList5<-list(as.data.frame(datasetsList3[[1]]),as.data.frame(datasetsList3[[2]]),as.data.frame(datasetsList3[[3]]),as.data.frame(datasetsListTopRank[[1]]),as.data.frame(datasetsListTopRank[[2]]),as.data.frame(datasetsListTopRank[[3]]),as.data.frame(datasetsListTopRank[[4]]))

m<-as.data.frame(datasetsListMeanRank[[1]])
n<-as.data.frame(datasetsListMultipliedPVal[[1]])
o<-as.data.frame(datasetsListMeanRank[[4]])

datasetsListMultipliedPValSmall<-combn(datasetsListSmall3,2,wrapperMethods,method=1)

write.table(m[,1],"~/forPresent4.tsv",col.names=T,row.names = T)

write.table(o,"~/datasetsListTopRank-Three.tsv",col.names = T)

## meanRankThree
meanRankThreeWrapper<-function(dataset1,dataset2,dataset3) {
  
  rows1<-gsub("-.*","",rownames(as.data.frame(dataset1)))
  rows2<-gsub("-.*","",rownames(as.data.frame(dataset2)))
  rows3<-gsub("-.*","",rownames(as.data.frame(dataset3)))
  
  together<-as.vector(sort(unique(union(rows1,union(rows2,rows3)))))
  
  test<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset1))),function(d) sapply(1:length(together),meanRankThree,b=d,dataset1=dataset1,dataset2=dataset2,dataset3=dataset3,together=together))))) #length(colnames(as.data.frame(dataset1))) #length(rownames(as.data.frame(newDataset)))
  colnames(test)<-colnames(dataset1)
  rownames(test)<-together
  return(test)
}

meanRankThree <-function(a,b,dataset1,dataset2,dataset3,together=together) {
  
    c<-paste0("^",together[a],"$")
    first<-grep(c,rownames(as.data.frame(dataset1)))
    second<-grep(c,rownames(as.data.frame(dataset2)))
    third<-grep(c,rownames(as.data.frame(dataset3)))
    
    if(!identical(first,integer(0))&&(!identical(second,integer(0))&&!identical(third,integer(0)))) {
      return(mean(c(as.data.frame(dataset1)[first,b],as.data.frame(dataset2)[second,b],as.data.frame(dataset3)[third,b])))
    }
    if(!identical(first,integer(0))&&(!identical(second,integer(0))&&identical(third,integer(0)))) {
      return(mean(c(as.data.frame(dataset1)[first,b],as.data.frame(dataset2)[second,b])))
    }
    if(!identical(first,integer(0))&&(identical(second,integer(0))&&!identical(third,integer(0)))) {
      return(mean(c(as.data.frame(dataset1)[first,b],as.data.frame(dataset3)[third,b])))
    }
    if(identical(first,integer(0))&&(!identical(second,integer(0))&&!identical(third,integer(0)))) {
      return(mean(c(as.data.frame(dataset2)[second,b],as.data.frame(dataset3)[third,b])))
    }
    if(!identical(first,integer(0))&&(identical(second,integer(0))&&identical(third,integer(0)))) {
      return((as.data.frame(dataset1))[first,b])
    } 
    if(identical(first,integer(0))&&(!identical(second,integer(0))&&identical(third,integer(0)))) {
      return((as.data.frame(dataset2))[second,b])
    } 
    if(identical(first,integer(0))&&(identical(second,integer(0))&&!identical(third,integer(0)))) {
      return((as.data.frame(dataset3))[third,b])
    } 
}
## Make intersection
datasetsListCombine<-list(chea2,encode2,coexp)

findMaxIntersect<-function(dataset1,dataset2) {
  intersection<-intersect(colnames(dataset1),colnames(dataset2))
  intersection<-intersection[!is.na(intersection)]
  maximum=0
  for(i in 1:length(intersection)) {
    a<-which(colnames(dataset1)==intersection[i])
    b<-which(colnames(dataset2)==intersection[i])
    print(a)
    summation1=character()
    summation2=character()
    for(j in 1:length(a)) {
      tf1<-as.data.frame(dataset1)[,a[j]]
      tf1<-tf1[!is.na(tf1)]
      tf1<-as.vector(tf1)
      summation1<-c(summation1,tf1)
    }
    for(k in 1:length(b)) {
      tf2<-dataset2[,b[k]]
      tf2<-tf2[!is.na(tf2)]
      tf2<-as.vector(tf2)
      summation2<-c(summation2,tf2)
    }
    summation<-unique(intersect(summation1,summation2))
    if(length(summation)>maximum) {
      maximum=length(summation)
    }
  }
  return(maximum)
}

wrapperIntersect<-function(list) {
  dataset1<-as.data.frame(list[[1]])
  dataset2<-as.data.frame(list[[2]])
  maximum<-findMaxIntersect(dataset1,dataset2)
  toAdd<-make_intersect_df(dataset1,dataset2,maximum)
  datasetsListCombine2[[length(datasetsListCombine2)+1]]<-toAdd
  return(datasetsListCombine2)
}

make_intersect_df<-function(dataset1,dataset2,maximum) {
  
  together<-sort(intersect(colnames(as.data.frame(dataset1)),colnames(as.data.frame(dataset2))))
  temp<-sapply(1:length(together),make_intersect,together=together,dataset1=dataset1,dataset2=dataset2,maximum=maximum)
  test<-cbind(temp[[1]],temp[[2]])
  for(i in 3:length(temp)){
    test<-cbind(test,temp[[i]])
  }
  test<-as.data.frame(test)
  # test<-as.data.frame(do.call(cbind,(lapply(1:length(together),make_intersect,together=together,dataset1=dataset1,dataset2=dataset2,maximum=maximum)))) #length(together)
  colnames(test)<-together
  return(test)
}

make_intersect<-function(b,together,dataset1,dataset2,maximum) {
  c<-paste0("^",together[b],"$")
  
  first<-grep(c,colnames(as.data.frame(dataset1)))
  second<-grep(c,colnames(as.data.frame(dataset2)))
  
  combined<-rep(NA,maximum)

  piece1<-character()
  piece2<-character()
  
  if(!identical(first,integer(0))) {
    for(i in 1:length(first)) {
      noNA<-as.data.frame(dataset1)[,first[i]]
      noNA<-noNA[!is.na(noNA)]
      noNA<-as.vector(noNA)
      piece1=c(piece1,noNA)
    }
    if(identical(second,integer(0))) {
      combined[1:length(piece1)]<-piece1
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
  }
  if(!identical(second,integer(0))) {
    for(j in 1:length(second)) {
      noNA<-as.data.frame(dataset2)[,second[j]]
      noNA<-noNA[!is.na(noNA)]
      noNA<-as.vector(noNA)
      piece2=c(piece2,noNA)
    }
    if(identical(first,integer(0))){
      combined[1:length(piece2)]<-piece2
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
  }
  if(!identical(first,integer(0))&&!identical(second,integer(0))) {
    piece<-unique(intersect(piece1,piece2))
 
    if(identical(piece,character(0))) {
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
    else {
      combined[1:length(piece)]<-piece
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
  }
}

temp<-list()
datasetsListCombine2<-list()
sample<-combn(datasetsListCombine,2,wrapperIntersect)
datasetsListCombine2<-sample

# m<-as.data.frame(datasetsList2[[1]])
# n<-as.data.frame(datasetsListCombine3[[1]])
# o<-as.data.frame(datasetsListCombine4[[2]])
# 
# write.table(o,"~/datasetsListCombine-three.tsv",col.names=T)

threeMax<-findMaxIntersect(as.data.frame(datasetsListCombine2[[1]]),as.data.frame(datasetsList[[3]]))
threeCombine<-make_intersect_df(as.data.frame(datasetsListCombine2[[1]]),as.data.frame(datasetsList[[3]]),threeMax)
datasetsListCombine2[[length(datasetsListCombine2)+1]]<-threeCombine #

## find p values
pValSimple<-function(tf1,tf2) {
  
  tf1<-tf1[!is.na(tf1)]
  tf2<-tf2[!is.na(tf2)]
  
  contTable<-matrix(c((20000-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  return(pVal[1])
}

makepVals<-function(x) {
  pVals<-as.data.frame(do.call(cbind,(apply(creedsFisher, 2, function(d) apply(x,2,pValSimple,tf1=d))))) #CHANGE TO CREEDSINTERSECTION
  
  pVals<-apply(pVals,c(1,2),unlist)
  pVals<-as.data.frame(pVals,stringsAsFactors=F)
  return(pVals)
}

ptm<-proc.time()
datasetsListCombine3<-lapply(datasetsListCombine2,makepVals)
proc.time()-ptm

datasetsListCombine4<-list(as.data.frame(datasetsList2[[1]]),as.data.frame(datasetsList2[[2]]),as.data.frame(datasetsList2[[3]]),as.data.frame(datasetsListCombine3[[1]]),as.data.frame(datasetsListCombine3[[2]]),as.data.frame(datasetsListCombine3[[3]]),as.data.frame(datasetsListCombine3[[4]]))
datasetsList5<-datasetsListCombine4
