## Reducing colnames to just tfs to make them easier to find
creeds_fisher_revised_colnames<-creedsFisher
colnames(creeds_fisher_revised_colnames)<-gsub("-.*","",colnames(creedsFisher))

new_creeds_revised_colnames<-new_creeds
colnames(new_creeds_revised_colnames)<-gsub("-.*","",colnames(new_creeds))

## Get all experiments with a certain tf from the datasets
datasets<-list(chea2,encode2,coexp,bioplex,humap)
CTCF<-find_exps("CTCF") 

JUN<-find_exps("JUN")

STAT3_old_creeds<-find_exps("STAT3",list(creeds_fisher_revised_colnames))
STAT3_exps<-find_exps("STAT3",datasets)

find_exps <-function(name,which_datasets) {
  exps<-list()
  for(i in 1:length(which_datasets)) {
    exps_dataset<-grep(paste0("^",name,"$"),colnames(as.data.frame(which_datasets[[i]])))
    # print(exps_dataset)
    if(!identical(exps_dataset,integer(0))) {
    for(j in 1:length(exps_dataset)) {
      with_nas<-which_datasets[[i]][,exps_dataset[j]]
      with_nas<-as.character(with_nas)
      exps[[length(exps)+1]]<-with_nas[!is.na(with_nas)]
    }
    }
  }
  return(exps)
}

## Using SuperExactTest package

jun_super<-supertest(JUN,n=20000)

plot(jun_super,degree = 2,sort.by = 'size',legend.col = 1)

plot(jun_super,Layout = 'landscape',degree = 2,sort.by = 'size') #11

write.csv(summary(jun_super)$Table, file="jun_super_creeds.csv", row.names=FALSE)

## Making table of intersection numbers

# Labeling the experiments. Can get dataset names and numbers of experiments from exp_distribution_df data frame (made in find_most_experiments.R)

l1<-make_label("chea",8) # will produce "chea1","chea2"..."chea8"
l2<-make_label("encode",3)

l1<-c(l1,l2)
l1[[length(l1)+1]]<-"bioplex1"

make_label <- function(dataset_name,num_exps) {
  l1<-list()
  for(i in 1:num_exps){
    l1[[length(l1)+1]]<-paste0(dataset_name,i)
  }
  return(l1)
}

## Get sizes of intersections between all of the datasets
a<-combn(JUN,2,function(x)return(length(Reduce(intersect,x))),simplify=T)

## Get sizes of intersections between creeds and all of the other datasets
stat3_df<-make_df(STAT3_old_creeds,STAT3_exps)

make_df <- function(creeds_exps,exps) {
  tf_df<-data.frame(matrix(NA, nrow=length(creeds_exps), ncol=length(exps)))
  for(i in 1:length(creeds_exps)) {
    row<-numeric()
    for(j in 1:length(exps)) {
      row<-c(row,length(intersect(creeds_exps[[i]],exps[[j]])))
    }
    tf_df[i,]<-row
  }
  rownames(tf_df)<-make_label("oldcreeds",length(creeds_exps)) 
  return(tf_df)
}

colnames(stat3_df)<-l1
