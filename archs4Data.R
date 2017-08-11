##set destination for raw expression data file
destination_file = "~/human_matrix_download.h5" # It might be human and mouse though.

# Check for dependencies and install if missing
packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
library("rhdf5")

# Check if gene expression file was already downloaded, if not in current directory download file from repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download.file(url, destination_file, quiet = FALSE)
} else{
  print("Local file already exists.")
}


## Retrieve information from compressed data ##
install.packages("bit64")

samples = h5read(destination_file, "meta/samples") # It said I couldn't open the file.
tissue = h5read(destination_file, "meta/tissue")
series = h5read(destination_file, "meta/series")
genes = h5read(destination_file, "meta/genes")
expression = h5read(destination_file, "data/expression")

## Quantile normalization
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')

library(preprocessCore)

#expression<-read.table("~/expression.tsv",sep="\t",stringsAsFactors = FALSE)
expression<-normalize.quantiles(expression)
write.table("~/normalizedCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)

pearsonCorr<-cor(expression)
write.table(pearsonCorr,"~/pearsonCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)

tfs<-read.delim("ADD/tfs.txt",header=FALSE,sep="\t")
tfsCoexp<-subset(pearsonCorr,colnames(pearsonCorr)%in%tfs)
write.table(tfsCoexp,"~/tfsCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)

absCoexp<-apply(tfsCoexp,c(1,2),abs(x)) # may have problems
write.table(absCoexp,"~/absCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)

cheaGenes<-read.delim("ADD/cheaGenes.txt",header=FALSE,sep="\t")
cheaCols<-read.delim("ADD/colnames(chea2).txt",header=FALSE,sep="\t")
cheaCoexp<-subset(absCoexp,colnames(absCoexp)%in%cheaCols)
cheaCoexp<-subset(cheaCoexp,rownames(cheaCoexp)%in%cheaGenes)
write.table(cheaCoexp,"~/cheaCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)

encodeGenes<-read.delim("ADD/encodeGenes.txt",header=FALSE,sep="\t")
encodeCols<-read.delim("ADD/colnames(encode2).txt",header=FALSE,sep="\t")
encodeCoexp<-subset(absCoexp,colnames(absCoexp)%in%encodeCols)
encodeCoexp<-subset(encodeCoexp,rownames(encodeCoexp)%in%encodeGenes)
write.table(encodeCoexp,"~/encodeCoexpression.tsv",sep="\t",row.names=TRUE,col.names=TRUE)
  
## get supplementary data
write(tfs,"~/tfs.txt",sep="\t")
write(colnames(chea2),"~/colnames(chea2).txt",sep="\t")
write(cheaGenes,"~/cheaGenes.txt",sep="\t")
write(colnames(encode2),"~/colnames(encode2).txt",sep="\t")
write(encodeGenes,"~/encodeGenes.txt",sep="\t")

cheaGenes<-unique(as.character(unlist(chea2)))
encodeGenes<-unique(as.character(unlist(encode2)))
