## to get a dataframe of how many experiments there are for each tf in each dataset.

bioplex<-read.table("~/bioplex.tsv",header = T,stringsAsFactors = F)
humap<-read.table("~/humap.tsv",header = T,stringsAsFactors = F)

tfs_colnames<-list(colnames(chea2),colnames(encode2),colnames(coexp),colnames(bioplex),colnames(humap),a,b)

a<-gsub("-.*","",colnames(creedsFisher))
b<-gsub("-.*","",colnames(new_creeds))

tfs<-Reduce(union,tfs_colnames)

exp_distribution_df<-most_exps(tfs)

exp_distribution_df[,2:9]<-apply(exp_distribution_df[,2:9],2,as.numeric)

write.table(exp_distribution_df,"exp_distribution_df.tsv",col.names = T)

most_exps <-function(tfs) {
  exp_distribution<-data.frame(tf=character(),chea=numeric(),encode=numeric(),coexp=numeric(),bioplex=numeric(),humap=numeric(),oldCreeds=numeric(),newCreeds=numeric(),total=numeric(),stringsAsFactors = F)
  for(i in 1:length(tfs)) { #
    summation = 0
    to_bind<-tfs[i]
    for(j in 1:length(tfs_colnames)) {
      
      c<-grep(paste0("^",tfs[i],"$"),tfs_colnames[[j]])
      if(j==6) {
       print(c) 
      }
      summation<-summation+length(c)
      to_bind<-c(to_bind,length(c))
    }
    to_bind<-c(to_bind,summation)
    to_bind<-as.data.frame(to_bind,stringsAsFactors = F)
    to_bind<-t(to_bind)
    colnames(to_bind)<-c("tf","chea","encode","coexp","bioplex","humap","oldCreeds","newCreeds","total")
    # print(to_bind)
    exp_distribution<-rbind(exp_distribution,to_bind,stringsAsFactors=F)
  }
  return(exp_distribution)
}
