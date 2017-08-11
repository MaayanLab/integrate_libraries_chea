## Import datasets, function from Ali
read.ttsv = function(file, header=TRUE, sep="\t", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

chea = read.ttsv("C:/Users/maayanlab1/Downloads/ChEA_2016.txt")
chea = chea[2:nrow(chea),]
encode = read.ttsv("C:/Users/maayanlab1/Downloads/ENCODE_TF_ChIP-seq_2015.txt")
encode = encode[2:nrow(encode),]

## Remove ", 1.0" from ChEA genes
chea1 = as.data.frame(gsub(",1.0","",as.matrix(chea)),stringsAsFactors = FALSE)

## Replace underscores in ChEA and ENCODE with dashes 
## (so column names are similar to those of creedsFisher and all can be benchmarked using the same code)
cheaDash<-chea1
colnames(cheaDash)<-gsub("_","-",colnames(chea1))

encodeDash<-encode
colnames(encodeDash)<-gsub("_","-",colnames(encode))

## Remove everything from column names except the TF name
chea2<-chea1
colnames(chea2)<-as.vector(gsub("_.*","",colnames(chea1)))

encode2<-encode
colnames(encode2)<-as.vector(gsub("_.*","",colnames(encode2)))

## Write to file
write.table(chea2,"C:/Users/maayanlab1/Downloads/chea2.tsv",sep="\t",col.names=TRUE)
write.table(encode2,"C:/Users/maayanlab1/Downloads/encode2.tsv",sep="\t",col.names=TRUE)

