# Need: chea2, encode2 

## Axes of the matrices. 
tfs<-union(colnames(chea2),colnames(encode2)) 
genes<-unique(as.character(union(unlist(chea3),unlist(encode2))))
genes = genes[!is.na(genes)]

# tfsBoth<-intersect(colnames(chea2),colnames(encode2)) 
# genesBoth<-intersect(unlist(chea3),unlist(encode2))
# genesBoth = genesBoth[!is.na(genesBoth)]


## Make AND adjacency matrix. 1 if gene is target of that tf, 0 if not. ChEA && ENCODE.
andNetwork<-matrix(0, length(genes),length(tfs))
andNetwork = as.data.frame(andNetwork)

colnames(andNetwork) = tfs
rownames(andNetwork) = genes

for(i in 1:length(tfs)) {
  chea_exp = chea2[,colnames(chea2)==tfs[i]]
  chea_exp = chea_exp[!is.na(chea_exp)]
  encode_exp = encode2[,colnames(encode2)==tfs[i]]
  encode_exp = encode_exp[!is.na(encode_exp)]
  for(j in 1:length(genes)){
    if(any(chea_exp == genes[j])&&any(encode_exp == genes[j])){
      andNetwork[j,i]=1
    }
  }
}

andNetwork1<-andNetwork[!is.na(row.names(andNetwork)),]
write.table(andNetwork,"C:/Users/maayanlab1/Downloads/andNetwork.csv",sep=",",row.names=TRUE,col.names=TRUE)

## Make OR adjacency matrix. 1 if gene is target of that tf in either ChEA or ENCODE, 0 if not. ChEA || ENCODE.
orNetwork<-matrix(0, length(genes),length(tfs))
orNetwork = as.data.frame(orNetwork)

colnames(orNetwork) = tfs
rownames(orNetwork) = genes

for(i in 1:length(tfs)) {
  chea_exp = chea2[,colnames(chea2)==tfs[i]]
  chea_exp = chea_exp[!is.na(chea_exp)]
  encode_exp = encode2[,colnames(encode2)==tfs[i]]
  encode_exp = encode_exp[!is.na(encode_exp)]
  for(j in 1:length(genes)){
    if(any(chea_exp == genes[j])||any(encode_exp == genes[j])){
      orNetwork[j,i]=1
    }
  }
}

write.table(orNetwork,"C:/Users/maayanlab1/Downloads/orNetwork.csv",sep=",",row.names=TRUE,col.names=TRUE)

## Calculating Jaccard distances
# install.packages("proxy")
library(proxy)

andDistances<-dist(andNetwork,method="Jaccard",by_rows=FALSE)
andDistances2<-as.matrix(andDistances)

orDistances<-dist(orNetwork,method="Jaccard",by_rows=FALSE)
orDistances2<-as.matrix(orDistances)

write.table(andDistances2,"C:/Users/maayanlab1/Downloads/andDistances2.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
write.table(orDistances2,"C:/Users/maayanlab1/Downloads/orDistances2.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
write.table(head(orDistances2),"C:/Users/maayanlab1/Downloads/headOrDistances2.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
d<-read.table("C:/Users/maayanlab1/Downloads/orDistances2.tsv",sep="\t")

## Jaccard, but only tfsBoth and genesBoth
tfsBoth<-intersect(colnames(chea2),colnames(encode2))
genesBoth<-intersect(unlist(chea3),unlist(encode2))
genesBoth = genesBoth[!is.na(genesBoth)]

andDistances3<-subset(andDistances2,colnames(andDistances2)%in%tfsBoth)
andDistances3<-subset(andDistances3,colnames(andDistances3)%in%tfsBoth)
