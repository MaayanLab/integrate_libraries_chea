## This methods is revised so that it can handle any number of datasets; not just 2.

datasetsListMultipliedPVal_ppi<-list()
datasetsListMeanRank_ppi<-list()
datasetsListTopRank_ppi<-list()

## executing each method
for(i in 2:length(datasetsList2_ppi)) {
  a<-combn(datasetsList2_ppi,i,wrapperMethods,method=1)
  for(j in 1:length(a)) {
    datasetsListMultipliedPVal_ppi[[length(datasetsListMultipliedPVal_ppi)+1]]<-a[[j]]
  }
}

ptm<-proc.time()
datasetsListMultipliedPVal_ppi_2<-combn(datasetsList2_ppi,2,wrapperMethods,method=1)
proc.time()-ptm

for(i in 2:length(datasetsList3_ppi)) {
  a<-combn(datasetsList3_ppi,i,wrapperMethods,method=2)
  for(j in 1:length(a)) {
    datasetsListMeanRank_ppi[[length(datasetsListMeanRank_ppi)+1]]<-a[[j]]
  }
}

for(i in 2:length(datasetsList3_ppi)) {
  a<-combn(datasetsList3_ppi,i,wrapperMethods,method=3)
  for(j in 1:length(a)) {
    datasetsListTopRank_ppi[[length(datasetsListTopRank_ppi)+1]]<-a[[j]]
  }
}


wrapperMethods<-function(list,method) {
  
  colsList<-list()
  for(i in 1:length(list)) {
    colsList[[length(colsList)+1]]<-gsub("-.*","",rownames(as.data.frame(list[[i]])))
  }
  together<-as.vector(sort(Reduce(union,colsList)))
  
  method_compute <-function(a,b,list,together,method) {
    c<-paste0("^",together[a],"$")
    grep_results<-list()
    
    for(i in 1:length(list)) {
      x<-grep(c,rownames(as.data.frame(list[[i]])))
      grep_results[[i]]<-as.data.frame(list[[i]])[x,b]
    }
    
    if(identical(method,1)){
      return(prod(unlist(grep_results),na.rm = T))
    }
    if(identical(method,2)){
      return(mean(unlist(grep_results),na.rm = T))
    }
    if(identical(method,3)){
      return(min(unlist(grep_results),na.rm = T))
    }
  }
  
  test<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(list[[1]]))),function(d) sapply(1:length(together),method_compute,b=d,list=list,together=together,method=method))))) 
  colnames(test)<-colnames(list[[1]])
  rownames(test)<-together
  
  datasetsList4[[length(datasetsList4)+1]]<-test
  return(datasetsList4)
}

