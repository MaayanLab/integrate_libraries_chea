## ENCODE. The col names without the GM
noGMEncode<-as.vector(gsub("_.*_","_",colnames(encode)))
encode3<-cbind(encode)
colnames(encode3)<-noGMEncode

## Find how many GMs were used
install.packages("gsubfn")
library(gsubfn)
middle<-as.vector(strapplyc(colnames(encode),"_(.*)_"))
GMs<-middle[grep('GM',middle)]
uniqueGMs<-unique(GMs)
length(unique(GMs))
length(unique(middle))

## Practice retrieving CREEDS IDs
sampleid<-character()
samplecount<-329
sampletf<-329
for(i in 1:2) {
  sampletf<-sampletf+1
  print(sampletf)
  
  r <- GET(paste0(CREEDS_URL, 'search'), query=list(q = tfs[i]))
  test <- fromJSON(httr::content(r, 'text'))
  if(length(test)!=0 && "hs_gene_symbol" %in% colnames(test)) {
    print("LOL")
    test = subset(test,test$hs_gene_symbol==tfs[i])
    test=subset(test,!(is.na(test["hs_gene_symbol"])))
    if(length(test$id)!=0) {
      sampleid<-append(sampleid, test$id)
      samplecount<-samplecount+1
    }
  }
  print("MEH")
}

## Make binary matrix, didn't reduce stuff to choose from so slower.
network<-matrix(0, length(genes),length(tfs))
network = as.data.frame(network)
colnames(network) = tfs
rownames(network) = genes
for(i in 1:length(tfs)) {
  chea_exp = chea3[,colnames(chea3)==tfs[i]]
  chea_exp = chea_exp[!is.na(chea_exp)]
  encode_exp = encode2[,colnames(encode2)==tfs[i]]
  encode_exp = encode_exp[!is.na(encode_exp)]
  for(j in 1:length(genes)){
    if(any(chea_exp == genes[j])&&any(encode_exp == genes[j])){
      network[j,i]=1
    }
    
  }
  
}

## Move to sources and targets table, first try
edges = data.frame("source","target")
colnames(edges) = c("x","y")
for(i in 1:length(colnames(network))){
  sources = rep(tfs[i],sum(network[,i]))
  targets = genes[network[,i]==1]
  temp_edges=data.frame(sources,targets)
  colnames(temp_edges)=c("x","y")
  edges = rbind(edges,temp_edges)
}

## Make network of tf-tf interactions. Did the same thing- called "reduced".
tfTargets<-subset(edges1, edges1[,2]%in%tfs)
tfTfInteractions<-write.table(tfTargets,file='tfTfInteractions.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

## To make reducedCREEDSCombined
count2<-0
for(i in 1) {
  count2<-count2+1
  print(count2)
  up<-stri_split_fixed(creeds[i,3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(creeds[i,4],",")
  down<-sapply(down,toupper)
  oneTF<-tfsInCREEDS[tfsInCREEDS[,1]==creeds[i,2],] # limit to the one tf of the row in creeds, with many targets
  if(length(oneTF[,1]) > 0) {
    upTF<-oneTF[oneTF[,2]%in%up,] # reduce to only targets in up regulated. maybe add info now, rbind.
    #print("problem? before")
    if(length(upTF[,1]) > 0) {
      upTF[,3]<-"up"
    }
    #print("problem after?")
    downTF<-oneTF[oneTF[,2]%in%down,]
    if(length(downTF[,1]) > 0) {
      downTF[,3]<-"down"
    }
    #print("here?")
    both<-rbind(upTF,downTF)
    #print("or here")
    if(length(both[,1]) > 0) {
      both[,4:6]<-NA
      both<-both[,c(1,3,2,4,5,6)]
      both[,4]<-creeds[i,1]
      both[,5]<-creeds[i,5]
      both[,6]<-creeds[i,6]
      
      colnames(both) = c("a","b","c","d","e","f")
      reducedCREEDSCombined<-rbind(reducedCREEDSCombined,both)
    }
    
  }
}

## Also for reducedCREEDSCombined
# Of the 75 transcription factors, only take the ones that are in CREEDS
reduced1<-reduced
reduced1[] <- lapply(reduced,as.character)
creeds1<-creeds

tfsInCREEDS<-reduced[reduced[,1]%in% creeds[,2],] #the tfs found in all 3 datasets
reducedCREEDSCombined<-data.frame("source", "interaction","target","series","cell type","perturbation", stringsAsFactors = FALSE)
colnames(reducedCREEDSCombined) = c("a","b","c","d","e","f")

count2<-0
for(i in 1:length(rownames(creeds))) {
  count2<-count2+1
  print(count2)
  up<-stri_split_fixed(creeds[i,3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(creeds[i,4],",")
  down<-sapply(down,toupper)
  oneTF<-tfsInCREEDS[tfsInCREEDS[,1]==creeds[i,2],] # limit to the one tf of the row in creeds, with many targets
  if(length(oneTF[,1]) > 0) {
    upTF<-oneTF[oneTF[,2]%in%up,] # reduce to only targets in up regulated. maybe add info now, rbind.
    #print("problem? before")
    if(length(upTF[,1]) > 0) {
      upTF[,3]<-"up"
    }
    #print("problem after?")
    downTF<-oneTF[oneTF[,2]%in%down,]
    if(length(downTF[,1]) > 0) {
      downTF[,3]<-"down"
    }
    #print("here?")
    both<-rbind(upTF,downTF)
    #print("or here")
    if(length(both[,1]) > 0) {
      both[,4:6]<-NA
      both<-both[,c(1,3,2,4,5,6)]
      both[,4]<-creeds[i,1]
      both[,5]<-creeds[i,5]
      both[,6]<-creeds[i,6]
      
      colnames(both) = c("a","b","c","d","e","f")
      reducedCREEDSCombined<-rbind(reducedCREEDSCombined,both)
    }
    
  }
}

## Test specific tf in CREEDS
response <- GET(paste0(CREEDS_URL, 'api'), query=list(id ='gene:1137'))
if (response$status_code == 200){
  response <- fromJSON(httr::content(response, 'text'))
  print(response)
}

## Labeling cell types
tFsCREEDSCellType[2:9,"g"]<-"hematopoietic cells"
tFsCREEDSCellType[10:11,"g"]<-"immune"
tFsCREEDSCellType[12:13,"g"]<-"liver"
tFsCREEDSCellType[14:18,"g"]<-"alveolar epithelial"
tFsCREEDSCellType[19,"g"]<-"immune"
tFsCREEDSCellType[20,"g"]<-"liver"
tFsCREEDSCellType[21:22,"g"]<-"immune"
tFsCREEDSCellType[23,"g"]<-"embryonic stem cells"
tFsCREEDSCellType[24,"g"]<-"neuron"
tFsCREEDSCellType[25:26,"g"]<-"immune"

## Getting a table for all tfs in CREEDS, and adding a column for human or mouse (organism)
install.packages("dplyr")
install.packages(c("nycflights13", "Lahman")) ## for examples

#tfsActuallyInCREEDS<-creeds[creeds[,2]%in% tfs,]
#tFsCREEDSCellType<-data.frame("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", stringsAsFactors = FALSE)

library(dplyr)
library(stringi)

tFsCREEDSCellType<-reducedCREEDSCombined
tFsCREEDSCellType[,c("g","h","i")]<-NA
tFsCREEDSCellType[1,7]<-"cell type"
tFsCREEDSCellType[1,8]<-"organism"
tFsCREEDSCellType[1,9]<-"dataset"
tFsCREEDSCellType[-1,9]<-"CREEDS"
tFsCREEDSCellType[1,4]<-"experiment (GEO series or PMID)"
tFsCREEDSCellType[1,5]<-"cell line"
tFsCREEDSCellType<-tFsCREEDSCellType[,c(1,2,3,4,5,7,6,8,9)]
#tFsCREEDSCellType[,7][tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedHuman[,4]]<-"human" another way to do it
tFsCREEDSCellType<-transform(tFsCREEDSCellType, h = ifelse(tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))

reducedCREEDSCombinedMouse<-anti_join(reducedCREEDSCombined,reducedCREEDSCombinedHuman)
tFsCREEDSCellType<-transform(tFsCREEDSCellType, h = ifelse(tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

## Get the cell lines I still need to identify
a<-creedsSep[!(creedsSep[,5]%in%allCellTypes[,1]),]

## Trying to process data so I can do it automatically
cellTypes2Revised<-read.table("C:/Users/maayanlab1/Downloads/cellTypes2 - cellTypes2.tsv",sep="\t",stringsAsFactors = FALSE)
uCT3Revised<-read.table("C:/Users/maayanlab1/Downloads/uniqueCellTypes3 - uniqueCellTypes3.tsv",sep="\t",stringsAsFactors = FALSE)

uCT3Revised[uCT3Revised[,1]=="Lin- Rosa26rtTA cells",2]<-"blood"

fill<-subset(uCT3Revised,is.na(uCT3Revised[,2])) # Make sure Lin-rosa had blood added; don't want to deal with spaces
fill1<-fill
fill[,1]<-sub(" (.*)","",fill[,1])

uCT3Revised<-uCT3Revised[!is.na(uCT3Revised[,2]),]

fillRevised<-read.table("C:/Users/maayanlab1/Downloads/fill - fill.tsv",sep="\t",stringsAsFactors = FALSE)
fillRevised[,1]<-fill1[,1]

allCellTypes<-rbind(cellTypes2Revised,uCT3Revised,fillRevised)
allCellTypes[719,1]<-creedsSepChEA[4596,5]

write.table(allCellTypes,"C:/Users/maayanlab1/Downloads/allCellTypes.tsv",sep="\t",row.names = FALSE,col.names = FALSE)

## Integrating cell types
cellTypes2<-read.csv("C:/Users/maayanlab1/Downloads/unique_cell_types - unique_cell_types.csv",stringsAsFactors = FALSE)
#cellTypes2<-as.data.frame(cellTypes2,stringsAsFactors=FALSE)
noNumbers<-as.vector(sub(".*? ","",cellTypes2[,1]))
cellTypes2[,1]<-noNumbers
for(i in 33:35) {
  cellTypes2[i,1]<-paste(cellTypes2[i,1],cellTypes2[i,3],sep=", ")
}
cellTypes2[,3]<-NULL
#cellTypes2<-cellTypes2[-234,] # it now skips row 234. 
colnames(cellTypes2)=c("cell lines","cell types")

cellTypes[,5]<-sub("TH1 ", "TH1",cellTypes[,5],fixed=TRUE) # space problem

for(i in 1:length(rownames(cellTypes2))) {
  cellTypes[which(cellTypes[,5]==cellTypes2[i,1]),6]<-as.character(cellTypes2[i,2])
  print(i)
}
NAs<-cellTypes[is.na(cellTypes[,6]),]

write.table(cellTypes2, file='C:/Users/maayanlab1/Downloads/cellTypes.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
write.table(cellTypes,file="C:/Users/maayanlab1/Downloads/cellTypesIn.csv",quote=FALSE,sep=',',row.names=FALSE,col.names=FALSE)

## Labeling GM cell types ENCODE
cellTypesNoGM<-cellTypes
cellTypesNoGM[,5][cellTypesNoGM[,5]=="GM18505"]<-"human"
cellTypesNoGM[,6][cellTypesNoGM[,5]=="GM18505"]<-"hematopoietic cells"
uniqueCellTypes<-as.data.frame(unique(cellTypes$`cell line`))

cellTypes<-read.table("C:/Users/maayanlab1/Downloads/cellTypesNoGMTSV.tsv", sep="\t",stringsAsFactors = FALSE)

## Parsing headers for ChEA
cellTypes<-tFsCREEDSCellType
# cellTypesCREEDSChEA<-data.frame("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", stringsAsFactors = FALSE)
# colnames(cellTypesCREEDSChEA)=c("a","b","c","d","e","f","g","h")

## Doesn't work as function for some reason, but separate seems to be okay.
parseHeaders<- function(smaller, cellTypes, tfs) {
  for(i in 1:length(colnames(smaller))) { #length(colnames(smaller))
    header<-stri_split_fixed(colnames(smaller)[i],"_")
    #column<-smaller[!is.na(smaller[,i]) ,i]
    #column<-"bo"
    column<-smaller[smaller[,i]%in%tfs,i]
    if(length(column)>0) {
      headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
      headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
      colnames(headerRow)=c("a","b","c","d","e","f","g","h")
      headerNotOneItemList<-unlist(header)
      if(length(headerNotOneItemList)!=3) { ## specifically to distinguish between chea and encode
        headerRow[,1]<-headerNotOneItemList[1]
        headerRow[,3]<-as.character(column)
        headerRow[,4]<-headerNotOneItemList[2]
        headerRow[,5]<-paste(headerNotOneItemList[5:length(headerNotOneItemList)-1],sep=" ",collapse = " ")
        headerRow[,8]<-headerNotOneItemList[length(headerNotOneItemList)]
      }
      else {
        headerRow[,1]<-headerNotOneItemList[1]
        headerRow[,3]<-column
        headerRow[,5]<-headerNotOneItemList[2]
        headerRow[,8]<-headerNotOneItemList[3]
      }
      print(i)
      print(headerRow)
      cellTypes<-rbind(cellTypes,headerRow)
    }
    
  }
}

parseHeaders(chea2,cellTypes,tfs)
parseHeaders(encode,cellTypes,tfs)

for(i in 1:length(colnames(chea2))) { #length(colnames(smaller))
  header<-stri_split_fixed(colnames(chea2)[i],"_")
  #column<-smaller[!is.na(smaller[,i]) ,i]
  #column<-"bo"
  column<-chea2[chea2[,i]%in%tfs,i]
  if(length(column)>0) {
    headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
    headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
    colnames(headerRow)=c("a","b","c","d","e","f","g","h","i")
    headerNotOneItemList<-unlist(header)
    
    headerRow[,1]<-headerNotOneItemList[1]
    headerRow[,3]<-as.character(column)
    headerRow[,4]<-headerNotOneItemList[2]
    headerRow[,5]<-paste(headerNotOneItemList[5:length(headerNotOneItemList)-1],sep=" ",collapse = " ")
    headerRow[,8]<-headerNotOneItemList[length(headerNotOneItemList)]
    headerRow[,9]<-"ChEA"
    
    print(i)
    print(headerRow)
    cellTypes<-rbind(cellTypes,headerRow)
  }
  
}
for(i in 1:length(colnames(encode))) { #length(colnames(smaller))
  header<-stri_split_fixed(colnames(encode)[i],"_")
  #column<-smaller[!is.na(smaller[,i]) ,i]
  #column<-"bo"
  column<-encode[encode[,i]%in%tfs,i]
  if(length(column)>0) {
    headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
    headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
    colnames(headerRow)=c("a","b","c","d","e","f","g","h","i")
    headerNotOneItemList<-unlist(header)
    
    headerRow[,1]<-headerNotOneItemList[1]
    headerRow[,3]<-column
    headerRow[,5]<-headerNotOneItemList[2]
    headerRow[,8]<-headerNotOneItemList[3]
    headerRow[,9]<-"ENCODE"
    
    print(i)
    print(headerRow)
    cellTypes<-rbind(cellTypes,headerRow)
  }
  
}

## Writing to file
colnames(creedsSepChEA)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target")
colnames(creedsSepENCODE)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target")
write.table(creedsSepChEA,"C:/Users/maayanlab1/Downloads/creedsSepChEA.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
write.table(creedsSepENCODE,"C:/Users/maayanlab1/Downloads/creedsSepENCODE.tsv",sep="\t",row.names = FALSE,col.names=TRUE)

## Learning how to get cell lines from online
#d<-read_html("http://web.expasy.org/cellosaurus/CVCL_0023")
#tmp <- tempfile(fileext = ".xml")
#write_xml(d, "C:/Users/maayanlab1/Downloads/alveolar.xml", options = "format")
#readLines(tmp)
#e<-htmlTreeParse("http://web.expasy.org/cellosaurus/CVCL_0023", useInternal=TRUE)

e<-getURL("http://web.expasy.org/cellosaurus/CVCL_0023",ssl.verifypeer=FALSE)
write(e,"C:/Users/maayanlab1/Downloads/alveolar.txt")

## Trying to actually use it:

#fill[,2]<-as.character(fill[,2])

## Trying to figure out functions
array<-matrix(0,2,2)
double = function(number,array) {
  number<-number*2
  array[1,2]<-number
  return(array)
}
array<-double(2, array)

## Just seeing if I could've made yesterday easier

cellLines_URL<-'https://scicrunch.org/resources/Cell%20Lines/search?q='
e<-getURL("http://web.expasy.org/cellosaurus/CVCL_0023",ssl.verifypeer=FALSE)
write(e,"C:/Users/maayanlab1/Downloads/alveolar.txt")

fail<-getURL("https://scicrunch.org/resources/Cell%20Lines/search?q=pre-ipscs&l=pre-ipscs",ssl.verifypeer=FALSE)
write(fail,"C:/Users/maayanlab1/Downloads/failedSearch.txt")


#cellLine<-fill[i,1]
# get_disease = function(cellLine, fill) {
#   assign('fill',fill,envir=.GlobalEnv)
#   cell_URL<- paste0('https://scicrunch.org/resources/Cell%20Lines/search?q=',cellLine,'&l=', cellLine, '&facet[]=Organism:Homo%20sapiens&facet[]=Organism:Mus%20musculus', collapse=NULL)
#   b<-getURL(cell_URL,ssl.verifypeer=FALSE)
#   write(b,"C:/Users/maayanlab1/Downloads/searchTest.txt")
#   searchTestFile<-read_file("C:/Users/maayanlab1/Downloads/searchTest.txt")
#   #print("a?")
#   if(!grepl("We could not find",searchTestFile)) {
#     #print("b?")
#     disease<-strapplyc(b,"Disease:([^\"]*)")
#     #print("c?")
#     disease<-unlist(disease)
#     #print(disease[1])
#     #print(i)
#     print(grep(cellLine, fill[,1]))
#     fill[grep(cellLine, fill[,1]),2]<-disease[1]
#     #print(fill[cellLine,2])
#     #return(fill)
#   }
#   
# }

## Trying to get URL
#a <- GET(paste0(GEO_URL,'gse2527','\\)&retmax=10&usehistory=y', collapse=NULL))
#b<-xmlParse(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL))
#b<-xmlParse(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL))
#b_data<-xmlToDataFrame(b)

## Notes:
## Get # tfs and targets in all three (function)
# getCount = function(dataset) { same as dim()
#   result<-c(length(colnames(dataset)),length(rownames(dataset)))
#   return(result)
# }

# colnames and rownames (or dim) can give. 
# chea tfs: 645, chea targets: 6208
# encode tfs: 816, encode targets: 6921

# They wanted all tf targets. 
# creeds tfs: 575, creeds targets: 

# get_organism = function(data_set_ids) {
#   GEO_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=('
#   b<-getURL(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL),ssl.verifypeer=FALSE)
# }

# for(i in 1:5) {
#   fill<-get_disease(fill[i,1],fill)
#   print(i)
# }

searchTest<-getURL(cell_URL,ssl.verifypeer=FALSE)
write(searchTest,"C:/Users/maayanlab1/Downloads/searchTest.txt")
searchTestFile<-read_file("C:/Users/maayanlab1/Downloads/searchTest.txt")

## Prepare to get whether mouse or human
b<-getURL(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL),ssl.verifypeer=FALSE)

query_key<-strapplyc(b,"<QueryKey>(.*)</QueryKey>")
web_env<-strapplyc(b,"<WebEnv>(.*)</WebEnv>")

second_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key='

c<-getURL(paste0(second_URL, query_key,'&WebEnv=',web_env, collapse=NULL),ssl.verifypeer=FALSE)
organism<-strapplyc(c,"Organism:\t(.*?)\nType:")
organism<-unlist(organism)
organism<-organism[1]

## Fill in cell type for some of them
toFill<-subset(creedsSep,is.na(creedsSep[,6])) # found out there are 345 cell types to fill.

creedsSepENCODE<-creedsSep[creedsSep[,3]%in%encodeTfs,]
toFill2<-subset(creedsSepENCODE,is.na(creedsSepENCODE[,6]))
uniqueCellTypes3<-unique(c(toFill2[,5],toFill1[,5]))
uniqueCellTypes3<-as.data.frame(uniqueCellTypes3)

uniqueCellTypes3[,2]<-NA

## Matching cell lines to predict cell types

for(i in 1:length(rownames(uniqueCellTypes3))) {
  if(any(grepl(paste(cellTypes2[,1],collapse="|"),uniqueCellTypes3[i,1],ignore.case=TRUE))) {
    print(cat("i=",i))
    for(j in 1:length(rownames(cellTypes2))) {
      print(cat("j=",j))
      if(grepl(paste(cellTypes2[j,1],collapse="|"),uniqueCellTypes3[i,1],ignore.case=TRUE)) {
        uniqueCellTypes3[i,2]<-cellTypes2[j,2]
        break
      }
    }
    #print(grep(paste(pattern,collapse="|"),uniqueCellTypes2[i,1],ignore.case=TRUE))
    #uniqueCellTypes2[i,2]<-cellTypes2[grep(paste(pattern,collapse="|"),uniqueCellTypes2[i,1],ignore.case=TRUE),2]
  }
}

##  Writing to file reducedCREEDSCombined
reducedCREEDSCombinedTSV<-write.table(reducedCREEDSCombined,file='C:/Users/maayanlab1/Downloads/reducedCREEDSCombinedTSV.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

## Removing non-human ones: 
reducedCREEDSCombinedHuman<-reducedCREEDSCombined[reducedCREEDSCombined[,4]!="GSE24594",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE2527",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE2433",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE40273",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE6846",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE38375",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE21060",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE1566",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE46970",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE16974",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE32224",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE39009",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE10954",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE4356",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE27159",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE55272",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE39443",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE30323",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE12999",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE31354",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE47989",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE18383",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE33659",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE45941",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE5654",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE19923",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE11664",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE1948",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE34545",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE40296",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE23923",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE28141",]

reducedCREEDSCombinedHumanTSV<-write.table(reducedCREEDSCombinedHuman,file='C:/Users/maayanlab1/Downloads/reducedCREEDSCombinedHumanTSV.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)



for(i in 1:length(rownames(creedsSepChEA))) {
  print(i)
  creedsSepChEA[i,6]<-allCellTypes[which(creedsSepChEA[i,5]==allCellTypes[,1]),2]
}


for(i in 1:length(rownames(creedsSepENCODE))) {
  creedsSepENCODE[i,6]<-as.character(allCellTypes[which(creedsSepChEA[i,5]==allCellTypes[,1]),2])
}

creedsSepChEA<-transform(creedsSepChEA, h = ifelse(creedsSepChEA[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))
creedsSepChEA<-transform(creedsSepChEA, h = ifelse(creedsSepChEA[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

creedsSepENCODE<-transform(creedsSepENCODE, h = ifelse(creedsSepENCODE[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))
creedsSepENCODE<-transform(creedsSepENCODE, h = ifelse(creedsSepENCODE[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

## Getting unique series from creedsSepChEA and creedsSepENCODE
noOrganismChEA<-subset(creedsSepChEA,is.na(creedsSepChEA[,8]))
noOrganismENCODE<-subset(creedsSepENCODE,is.na(creedsSepENCODE[,8]))
noOrganism<-unique(rbind(noOrganismChEA,noOrganismENCODE))

## creedsSep after \t halted it
for(i in 667516:length(rownames(creedsSep))) {
  print(i)
  creedsSep[i,6]<-allCellTypes[which(creedsSep[i,5]==allCellTypes[,1]),2] 
}

## New matrix. The reduced stuff.
tfsBoth<-intersect(colnames(chea3),colnames(encode2))
genesBoth<-as.character(intersect(unlist(chea3),unlist(encode2)))
genesBoth = genesBoth[!is.na(genesBoth)]

reducedNetwork<-matrix(0, length(genesBoth),length(tfsBoth))
reducedNetwork = as.data.frame(reducedNetwork)
colnames(reducedNetwork) = tfsBoth
rownames(reducedNetwork) = genesBoth
for(i in 1:length(tfsBoth)) {
  chea_exp = chea3[,colnames(chea3)==tfsBoth[i]]
  chea_exp = chea_exp[!is.na(chea_exp)]
  encode_exp = encode2[,colnames(encode2)==tfsBoth[i]]
  encode_exp = encode_exp[!is.na(encode_exp)]
  for(j in 1:length(genesBoth)){
    if(any(chea_exp == genes[j])&&any(encode_exp == genes[j])){
      reducedNetwork[j,i]=1
    }
  }
}
# So in reducedNetwork, all tfs are in both datasets and all targets are too. But targets are not only tfs.
# It's reduced that has targets that are just tfs.

# Only have targets (genes) that are tfs
reducedEdges = data.frame("source","target")
colnames(reducedEdges) = c("x","y")
for(i in 1:length(colnames(reducedNetwork))){
  sources = rep(tfsBoth[i],sum(reducedNetwork[,i]))
  targets = genesBoth[reducedNetwork[,i]==1]
  temp_reducedEdges=data.frame(sources,targets)
  colnames(temp_reducedEdges)=c("x","y")
  reducedEdges = rbind(reducedEdges,temp_reducedEdges)
}

reduced<-subset(reducedEdges, reducedEdges[,2]%in%tfs) # see above note. Reduced was made with edges in ChEA or ENCODE though. Not &&
reducedTSV<-write.table(reduced,file='C:/Users/maayanlab1/Downloads/reducedTSV.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

## Example of calculating Jaccard distance
a<-c(0,1,1)
b<-c(0,1,0)
matrix<-cbind(a,b)
dist(matrix,method="Jaccard",by_rows=FALSE)
a[2]<-1

## Even simpler function that still doesn't work

b = function(a,c) {
  return(a+c)
}

c<-b(2,5)

## Before I turned it into a function.
tfsCreedsChEA<-unique(creedsSepChEAHuman[,1])
for(i in 1:length(tfsCreedsChEA)) {
  act<-sum(creedsSepChEAHuman[,1]==tfsCreedsChEA[i]&((creedsSepChEAHuman[,7]%in%activation&creedsSepChEAHuman[,2]=="up")|(creedsSepChEAHuman[,7]%in%inhibition&creedsSepChEAHuman[,2]=="down")))
  inhib<-sum(creedsSepChEAHuman[,1]==tfsCreedsChEA[i]&((creedsSepChEAHuman[,7]%in%activation&creedsSepChEAHuman[,2]=="down")|(creedsSepChEAHuman[,7]%in%inhibition&creedsSepChEAHuman[,2]=="up")))
  creedsSepChEAHuman[creedsSepChEAHuman[,1]==tfsCreedsChEA[i],10]<-act
  creedsSepChEAHuman[creedsSepChEAHuman[,1]==tfsCreedsChEA[i],11]<-inhib
  
  up<-data.frame(tfsCreedsChEA[i],act,inhib,"creedsSepChEAHuman",stringsAsFactors = FALSE)
  colnames(up)=c("a","b","c","d")
  upDownCounts<-rbind(upDownCounts,up)
  
}

tfsCreedsENCODE<-unique(creedsSepENCODEHuman[,1])
for(i in 1:length(tfsCreedsENCODE)) {
  act<-sum(creedsSepENCODEHuman[,1]==tfsCreedsENCODE[i]&((creedsSepENCODEHuman[,7]%in%activation&creedsSepENCODEHuman[,2]=="up")|(creedsSepENCODEHuman[,7]%in%inhibition&creedsSepENCODEHuman[,2]=="down")))
  inhib<-sum(creedsSepENCODEHuman[,1]==tfsCreedsENCODE[i]&((creedsSepENCODEHuman[,7]%in%activation&creedsSepENCODEHuman[,2]=="down")|(creedsSepENCODEHuman[,7]%in%inhibition&creedsSepENCODEHuman[,2]=="up")))
  creedsSepENCODEHuman[creedsSepENCODEHuman[,1]==tfsCreedsENCODE[i],10]<-act
  creedsSepENCODEHuman[creedsSepENCODEHuman[,1]==tfsCreedsENCODE[i],11]<-inhib
  
  up<-data.frame(tfsCreedsENCODE[i],act,inhib,"creedsSepENCODEHuman",stringsAsFactors = FALSE)
  colnames(up)=c("a","b","c","d")
  upDownCounts<-rbind(upDownCounts,up)
}

## Try
uniqueTfs<-unique(creedsSepChEAHuman[,1])
for(i in 1:length(uniqueTfs)) {
  creedsSepChEAHuman[(creedsSepChEAHuman[,1]==uniqueTfs[i]&((creedsSepChEAHuman[,7]%in%activation&creedsSepChEAHuman[,2]=="up")|(creedsSepChEAHuman[,7]%in%inhibition&creedsSepChEAHuman[,2]=="down"))),2]<-"activation"
  print("A")
  creedsSepChEAHuman[(creedsSepChEAHuman[,1]==uniqueTfs[i]&((creedsSepChEAHuman[,7]%in%activation&creedsSepChEAHuman[,2]=="down")|(creedsSepChEAHuman[,7]%in%inhibition&creedsSepChEAHuman[,2]=="up"))),2]<-"repression"
  print("b")
}

# updated<-list.files(path="~/sig_files_7-6-17",pattern=".*_") 
# notUpdated<-setdiff(files,updated)

## sep function. Sad I couldn't use it.
sep<-function(file,ppi) {
  
  if(grepl("2017",basename(file))) {
    print(basename(file))
    extract<-read.table(paste0("~/sig_files_7-6-17/",file),sep="\t",stringsAsFactors = FALSE)
    #print(head(extract))
    #print("break")
    ppi<-extractPPI(extract,ppi,basename(file))
    #assign("ppi",ppi,.GlobalEnv)
    # print("ppi in if statement")
    # print(head(ppi))
    # print("end of if")
    return(ppi)
  }
  else {
    print(basename(file))
    extract<-read.table(paste0("~/sig_files_7-6-17/",file),sep=" ",stringsAsFactors = FALSE)
    ppi<-extractPPI(extract,ppi,basename(file))
    #assign("ppi",ppi,.GlobalEnv)
    # print("ppi in else statement")
    # print(head(ppi))
    # print("end of else \n")
    return(ppi)
  }
  
}

ppi<-sapply(files,sep,ppi=ppi)
ppi1<-ppi
as.data.frame(ppi)

## cheaTfs, learning piping
cheaTfs<-unlist(chea2)
cheaTfs<-subset(cheaTfs,cheaTfs%in%tfs) #Just the targets in chea that are tfs.
cheaTfs<-union(colnames(chea2),cheaTfs) # anything in chea that is tf: colname and in data frame ## wanted UNION; would have had to use unique(c(,))

creedsSepENCODE<-creedsSep[creedsSep[,1]%in%colnames(encode2),]
encodeTfs<-unlist(encode2)
encodeTfs<-subset(encodeTfs,encodeTfs%in%tfs) 
encodeTfs<-union(colnames(encode2),encodeTfs)

creedsSepENCODE<-creedsSepENCODE[creedsSepENCODE[,3]%in%encodeTfs,]

creedsSepChEA<-creedsSep[creedsSep[,1]%in%colnames(chea2),]
columnNames <- creedsSep %>%
  rename(.,c("a banana"=V1,"2 bananas"=V2))
rename(creedsSep, V1='a')

## Narrow to targets that are also in chea and are tfs
creedsSepChEA<-creedsSepChEA[creedsSepChEA[,3]%in%cheaTfs,]
creedsSepENCODE<-creedsSepENCODE[creedsSepENCODE[,3]%in%encodeTfs,]

## Make only human
creedsSepChEAHuman<-creedsSepChEA[creedsSepChEA[,8]=="Homo sapiens",]
creedsSepENCODEHuman<-creedsSepENCODE[creedsSepENCODE[,8]=="Homo sapiens",] 

## Export
colnames(creedsSepChEAHuman)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target","activations","repressions")
colnames(creedsSepENCODEHuman)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target","activations","repressions")

write.table(creedsSepChEAHuman,"~/creedsSepChEAHuman.tsv",sep="\t",row.names=FALSE, col.names=TRUE)
write.table(creedsSepENCODEHuman,"~/creedsSepENCODEHuman.tsv",sep="\t",row.names=FALSE, col.names=TRUE)

cheaTfs %<>% unlist(chea2) %>% subset(.,.%in%tfs) %>% union(colnames(chea2),.)

creedsSepChEA1<- creedsSep %>% 
  filter(.,V1%in%colnames(chea2)) %>% # Narrow to sources are also tfs in chea
  filter(.,V3%in%cheaTfs) %>% # Narrow to targets that are also in chea and are tfs
  filter(.,V8 == "Homo sapiens") %>% # Make only human
  rename(.,source=V1,interaction=V2,target=V3,"experiment (GEO series or PMID)"=V4,"cell line"=V5, "cell type"=V6, "perturbation"=V7, "organism"=V8, "source:target"=V9) %T>%
  write.table(.,"~/pipe_test.tsv",sep="\t",row.names=FALSE, col.names=TRUE)

## Same for ENCODE: 
encodeTfs %<>% unlist(encode2) %>% subset(.,.%in%tfs) %>% union(colnames(encode2),.)

#colnames(upDownCounts)=c("transcription factor","# activations","# repressions","dataset")
upDownCounts<-upDownCounts[-1,]

write.table(creedsTfs,"C:/Users/maayanlab1/Downloads/creedsTfs.tsv",sep="\t",row.names = FALSE,col.names=TRUE)

## Try
a<-matrix(c(1,2,3,4),2,2)
colnames(a)=c("a","b")
b<-subset(a,colnames(a)=="a")
c<-subset(a,colnames(a)=="a")

# repeats<-ppiTfs[,1:2]
# ppiTfs0<-ppiTfs[!duplicated(repeats),]
# ppiTfs0<-ppiTfs[(duplicated(ppiTfs[,1:2]) | duplicated(ppiTfs[,1:2], fromLast = TRUE)), ]
# ppiTfsDiff<-apply(ppiTfs,1,!(%in%ppiTfs0))

fruits<-read.table("~/fruits.tsv",sep="\t")
colnames(fruits)<-unlist(fruits[1,])
fruits<-fruits[-1,]
fruits <- data.frame(matrix(unlist(fruits), nrow=5, byrow=T),stringsAsFactors=FALSE) # not exactly what I wanted but oh well
fruitsTable<-matrix(0,8,3)
fruitsTable<-as.data.frame(fruitsTable)

basket<-c("apple","cherry","sprout","carrot","elderberry")
eat<-NULL
fruitsTable<-apply(fruits,2,function(i) drink(i))

drink<- function(i) {
  print(i)
  targets<-match(i,basket)
  print(targets)
  targets<-targets[!is.na(targets)]
  print(targets)
  fruitsTable[targets,i]<-1
}

## These likely won't work because when you use apply it takes the column as a vector. It is viewed separately from the data frame.
cheaEditedNetwork1<-sapply(1:length(b),function(i) cheaEditedNetwork1[b[i],1]<-1)
cheaEdited1<-apply(cheaEdited,2,function(x){x[is.na(x)] <- 0 ; x})
cheaEditedNetwork2<-apply(cheaEdited,2,function(i) findColumn,cheaEdited=cheaEdited,cheaEditedNetwork1,cheaEditedNetwork1)

#cheaEditedNetwork[cheaEditedNetwork[,1]%in%b,1]<-1

## Nice but not actually what we wanted
cheaPVal<-matrix(data=NA,nrow=length(colnames(chea)),ncol=length(colnames(chea)))
rownames(cheaPVal)<-colnames(chea)
colnames(cheaPVal)<-colnames(chea)

## For 2 datasets
contTable<-matrix(c((length(union(dataset1,dataset2))-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(col1,col2))),nrow=2)

## Figuring out getting rid of NAs
a<-chea1[["FOXP1_22492998_ChIP.Seq_STRATIUM_Mouse"]]
a<-as.vector(a)
a[!is.na(a)]

# a<-dataset[[d[1]]]
# a<-as.vector(a)
# a<-a[!is.na(a)]
# print(a)
# print(typeof(table))

## figuring out r bind with function. I guess apply produces something each time, which you bind outside of the function.
example<-function(row,bind){
  # row<-data.frame("a","b","c",stringsAsFactors = FALSE)
  colnames(row)=c("1","2","3")
  rbind(bind,row)
}

bind=data.frame(apple=character(),banana=character(),cherry=character(),stringsAsFactors = FALSE)

bind<-lapply(l,example,bind=bind)

fishers<-function(d,dataset,table) {
   print(d)
   #print("Hello")
   tf1<-dataset[[d[1]]]
   print(dataset[[d[1]]])
   tf1<-tf1[!is.na(tf1)]

   tf2<-dataset[[d[2]]]
   tf2<-tf2[!is.na(tf2)]

   contTable<-matrix(c((length(dataset)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
   pVal<-fisher.test(contTable)
   table<-rbind(table,c(d[1],d[2],pVal))
   return(table)
 }

 cheaPVal<-lapply(c,function(d) fishers(d),dataset=chea1,table=cheaPVal)

d<-chea1["FOXP1_22492998_ChIP.Seq_STRATIUM_Mouse"]
d<-chea1$FOXP1_22492998_ChIP.Seq_STRATIUM_Mouse
d[!is.na(d)]

sapply(test[1:10,1],targetGenes) #length(chea[x])

## dplyr solution
library(dplyr)
ptm <- proc.time()
cheaPVal<-lapply(c, fishers,dataset=chea1) %>% bind_rows()

## function to test stuff on
example<-function(row,bind){
  # row<-data.frame("a","b","c",stringsAsFactors = FALSE)
  colnames(row)=c("1","2","3")
  rbind(bind,row)
}

bind=data.frame(apple=character(),banana=character(),cherry=character(),stringsAsFactors = FALSE)

bind<-lapply(l,example,bind=bind)

## Trying to make table with the 4 contingency values

test$inAinB<-NA
test$notBinA<-NA
test$notAinB<-NA
test$notAnotB<-NA

targetGenes<-function(x,genes,test) {
  #print(x)
  #print("a")
  print(test[x,1])
  print("hey")
  col1<-chea1[test[x,1]]
  
  col1<-col1[!is.na(col1)]
  #print(col1)
  #print("c")
  
  col2<-chea1[test[x,2]]
  col2<-col2[!is.na(col2)]
  #print(col2)
  
  #intersect<-length(intersect(col1,col2))
  test[x,3]<-length(intersect(col1,col2))
  
  # AB1<-length(col1)-intersect
  # AB2<-length(col2)-intersect
  test[x,4]<-length(col1)-test[x,3]
  test[x,5]<-length(col2)-test[x,3]
  
  test[x,6]<-length(genes)+test[x,3]-test[x,4]-test[x,5]
  return(test)
  # print(x[1])
  # print("pause")
  # print(x[2])
  # print("end")
  # return(length(col))
  #print(chea[x])
}

mini<-test[1:2,]
mini<-as.data.frame(mini,stringsAsFactors = FALSE)
#test$tf1Length<-NA
mini<-sapply(1:length(rownames(mini)),targetGenes,genes=genes,test=mini) #test$intersect[1:5]

test$tf2Length<-sapply(test$exp2,targetGenes)

col<-chea1["BP1_19119308_ChIP.ChIP_Hs578T_Human"]
col<-col[!is.na(col)]

col<-as.vector(col)
length(col)

b<-test[c("exp1","exp2")]

test$tf2Length<-NA
test$pVal<-NA

#cheaColsPrev<-combn(colnames(chea1),2,simplify=FALSE)
test<-as.data.frame(cheaColsComboArray,stringsAsFactors = FALSE) 
colnames(test)=c("exp1","exp2")

## Fisher p value stuff as applied to ChEA
cheaColsPrev<-combn(colnames(chea1),2,simplify=FALSE)
#cheaColsCombo<-t(combn(colnames(chea1),2,simplify=TRUE))

c<-cheaColsPrev[1:100]

ptm <- proc.time()
cheaPVal<-data.frame(tf1=character(),tf2=character(),pVal=character(),stringsAsFactors = FALSE)
genes<-unique(unlist(chea1))
genes<-genes[!is.na(genes)]


fishers<-function(d) { # will have to change function for comparing chea and encode
  
  tf1<-chea1[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  #print(cat(tf1,"hello"))
  
  tf2<-chea1[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  #print(tf2)
  
  #print(d)
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  #print(contTable)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  #print(pVal[[1]])
  row<-data.frame(d[1],d[2],pVal[[1]],stringsAsFactors = FALSE)
  
  colnames(row)=c("tf1","tf2","pVal")
  return(row)
  
}

cheaPVal<-as.data.frame(do.call(rbind,(lapply(cheaColsPrev, function(d) fishers(d)))))

proc.time() - ptm

write.table(cheaPVal,"~/cheaPVal.tsv",col.names=TRUE,sep="\t")

## turn other way
encodeColsCombo<-combn(colnames(encode),2,simplify=FALSE)

encodeColsSimple<-combn(colnames(encode2),2,simplify=TRUE)
encodeColsCombo<-t(combn(colnames(encode),2,simplify=TRUE))

## example to test setdiff with lists
both<-combn(union(colnames(chea2),c("c","d")),2,simplify=FALSE)
# ab<-combn(,2,simplify=FALSE)
# cd<-combn(c("c","d"),2,simplify=FALSE)
both<-setdiff(both,ab)
both<-setdiff(d,cd)

  bothColsCombo<-combn(union(colnames(chea1),colnames(encode)),2,simplify=FALSE)
bothColsCombo<-setdiff(bothColsCombo,cheaColsPrev)

test<-t(combn(union(colnames(chea1),colnames(encode)),2,simplify=TRUE))
test<-as.data.frame(test,stringsAsFactors=FALSE)

encodeComboDF<-t(combn(colnames(encode),2,simplify=TRUE))
encodeComboDF<-as.data.frame(encode,stringsAsFactors=FALSE)

cheaComboDF<-t(combn(colnames(chea1),2,simplify=TRUE))
cheaComboDF<-as.data.frame(cheaComboDF,stringsAsFactors=FALSE)

#get rid of ones in those

bothColsCombo<-setdiff(bothColsCombo,union(encodeColsCombo,cheaColsPrev))

res <- lapply(a, function(ch) grep("CREB1_HepG2_hg19", ch))

# which vectors contain a search term
b<-sapply(res, function(x) length(x) > 0)

## to see:
a<-t(bothColsCombo)


a<-data.frame(exp1=character(526320),exp2=character(526320),stringsAsFactors = FALSE)
b<-lapply(bothColsCombo,function(x) { return(x[1]) })
b<-unlist(b)
a[,1]<-b

c<-lapply(bothColsCombo,function(x) { return(x[2]) })
c<-unlist(c)
a[,2]<-c

bothColsComboDF<-a
write.table(bothColsComboDF,"~/bothColsComboDF.tsv",col.names=T,sep="\t")

## used for bothColsCombo
bothColsCombo<-combn(union(colnames(chea1),colnames(encode)),2,simplify=FALSE)
bothColsCombo<-setdiff(bothColsCombo,cheaColsPrev)
bothColsCombo<-setdiff(bothColsCombo,encodeColsCombo)

c<-bothColsCombo[1:100]
d<-c

genes<-unique(union(unlist(encode),unlist(chea1)))
genes<-genes[!is.na(genes)]

##
ptm <- proc.time()

b<-lapply(bothColsCombo,function(x){
  
  if(any(grepl(x[1],colnames(encode)))){
    #print(x)
    
    #x<-rev(x)
    return(x)
  }
  #return(x)
  
})
proc.time() - ptm
##
bothPVal<-data.frame(tf1=character(),tf2=character(),pVal=character(),stringsAsFactors = FALSE)
ptm <- proc.time()

fishers<-function(d) { # will have to change function for comparing chea and encode
  
  tf1<-chea1[[d[1]]]
  if(length(tf1)==0){
    print(d)
  }
  tf1<-tf1[!is.na(tf1)]
  #print(cat(tf1,"hello"))
  
  tf2<-encode[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  #print(tf2)
  
  #print(d)
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  #print(contTable)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  #print(pVal[[1]])
  row<-data.frame(d[1],d[2],pVal[[1]],stringsAsFactors = FALSE)
  
  colnames(row)=c("tf1","tf2","pVal")
  return(row)
  
}

bothPVal<-as.data.frame(do.call(rbind,(lapply(bothColsCombo, function(d) fishers(d)))))

proc.time() - ptm

write.table(bothPVal,"~/bothPVal.tsv",col.names=TRUE,sep="\t")
#encodePVal<-read.table("~/encodePVal.tsv",header=TRUE,stringsAsFactors = FALSE)

## The code I used to make CREEDS. Now revised to get rid of expression values.
CREEDS_URL <- 'http://amp.pharm.mssm.edu/CREEDS/'

ids <- character()

  r <- GET(paste0(CREEDS_URL, 'search'), query=list(q = "FOXP1"))
  response <- fromJSON(httr::content(r, 'text'))
  if(length(response)!=0 && "hs_gene_symbol" %in% colnames(response)) {
    response = subset(response,response$hs_gene_symbol==tfs[i])
    response<-subset(response,!(is.na(response["hs_gene_symbol"])))
    if(length(response$id)!=0) {
      ids<-append(ids, response$id)
    }
  }


## Make table
creeds<-matrix(NA,count,6)
creeds<-as.data.frame(creeds)

count1=0
## Put into table
for(i in 1:length(ids)) {
  count1<-count1+1
  print(count1)
  r <- GET(paste0(CREEDS_URL, 'api'), query=list(id= 1809))
  sig <- fromJSON(httr::content(r, 'text'))
  creeds[i,1]<-sig$geo_id
  creeds[i,2]<-sig$hs_gene_symbol
  creeds[i,3]<-paste0(sig$up_genes,collapse=",")
  creeds[i,4]<-paste0(sig$down_genes,collapse=",")
  creeds[i,5]<-sig$cell_type
  creeds[i,6]<-sig$pert_type
  
}

##
count1=0
## Put into table
for(i in 1:length(ids)) {
  count1<-count1+1
  print(count1)
  r <- GET(paste0(CREEDS_URL, 'api'), query=list(id= ids[i]))
  sig <- fromJSON(httr::content(r, 'text'))
  creeds[i,1]<-sig$geo_id
  creeds[i,2]<-sig$hs_gene_symbol
  creeds[i,3]<-paste0(sig$up_genes,collapse=",")
  creeds[i,4]<-paste0(sig$down_genes,collapse=",")
  creeds[i,5]<-sig$cell_type
  creeds[i,6]<-sig$pert_type
  
}

## Revised table:
# creeds1<-matrix(NA,count,6) didn't work because forgot to change creeds[i,1]<-sig$geo_id
# creeds1<-as.data.frame(creeds1)

## Old way to make creedsSep
count2<-0
for(i in 1:length(rownames(creeds))) { #:length(rownames(creeds))
  count2<-count2+1
  print(count2)
  up<-stri_split_fixed(creeds[i,3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(creeds[i,4],",")
  down<-sapply(down,toupper)
  
  sourcesUp = rep(creeds[i,2],length(up))
  interactionsUp = rep("up",length(up))
  tempUp<-data.frame(sourcesUp,interactionsUp,up)
  colnames(tempUp) = c("a","b","c")
  
  sourcesDown = rep(creeds[i,2],length(down))
  interactionsDown = rep("down",length(down))
  tempDown<-data.frame(sourcesDown,interactionsDown,down)
  colnames(tempDown) = c("a","b","c")
  
  creedsSepTemp = rbind(tempUp,tempDown)
  
  creedsSepTemp[,4:9]<-NA
  creedsSepTemp[,4]<-creeds[i,1]
  creedsSepTemp[,5]<-creeds[i,5]
  #creedsSepTemp[,6]<-cellTypes2[which(creedsSepTemp[i,5]==cellTypes2[,1]),2]
  creedsSepTemp[,7]<-creeds[i,6]
  creedsSepTemp[i,9]<-paste0(creedsSepTemp[i,1],":",creedsSepTemp[i,3],collapse=NULL) # may see problems here # forgot [i,9]-should be fine now.
  colnames(creedsSepTemp) = c("a","b","c","d","e","f","g","h","i")
  creedsSep<-rbind(creedsSep,creedsSepTemp)
  
}

## putting cellTypes in, the old for loop way
for(i in 1:length(rownames(creedsSep))) { 
  #print(i)
  creedsSep[i,6]<-allCellTypes[which(creedsSep[i,5]==allCellTypes[,1]),2]
}

## fix source:target, the old for loop way
for(i in 1:length(rownames(creedsSep))) {
  creedsSep[i,9]<-paste0(creedsSep[i,1],":",creedsSep[i,3],collapse=NULL) 
}

## Seems like for loop did better than apply
a<-apply(creedsSep[1:10000,],1,function(x){return(creedsSep1[which(x[9]==creedsSep1[,9])[1],8])})
proc.time() - ptm
creedsSep2<-creedsSep
creedsSep2[,9]<-a

genes<-unique(unlist(creedsFisher))
genes<-genes[!is.na(genes)]

genes<-unique(union(unlist(chea1),unlist(creedsFisher)))
genes<-genes[!is.na(genes)] #also time considerations though

if(dataset1==dataset2) {
  if(dataset==creedsFisher){
    write.table(datasetPVal,"~/creedsPVal.tsv",col.names=TRUE,sep="\t") 
    break
    if(dataset==chea1) {
      write.table(datasetPVal,"~/cheaPVal.tsv",col.names=TRUE,sep="\t")
    }
    else {
      write.table(datasetPVal,"~/encodePVal.tsv",col.names=TRUE,sep="\t")
    }
  }
}

## encode fisher
encodeColsCombo<-combn(colnames(encode),2,simplify=FALSE)
c<-encodeColsCombo[1:100]
genes<-unique(unlist(encode))%>%na.omit()

ptm <- proc.time()

fishersJaccard<-function(d) { 
  
  tf1<-encode[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  
  tf2<-encode[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  
  jaccard<-(1-(length(intersect(tf1,tf2))/length(union(tf1,tf2))))
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  row<-data.frame(d[1],d[2],pVal[[1]],jaccard,length(genes),length(tf1),length(tf2),length(intersect(tf1,tf2)),stringsAsFactors = FALSE)
  colnames(row)=c("tf1","tf2","pVal","jaccard","background genes","tf1 target genes","tf2 target genes","intersection size")
  return(row)
}

encodePVal<-as.data.frame(do.call(rbind,(lapply(encodeColsCombo, function(d) fishersJaccard(d))))) # PRACTICE

proc.time() - ptm

write.table(encodePVal,"~/encodePVal.tsv",col.names=TRUE,sep="\t")

encodePVal$match<-NA
encodePVal[["match"]]<-apply(encodePVal,1,function(x) {ifelse(gsub("_.*","",x[1])==gsub("_.*","",x[2]),1,0)})
write.table(encodePVal,"~/encodePValJMatch.tsv",col.names=TRUE,sep="\t")

## chea fisher
cheaColsCombo<-combn(colnames(chea1),2,simplify=FALSE)
c<-cheaColsCombo[1:100]
genes<-unique(unlist(chea1))%>%na.omit()

ptm <- proc.time()

fishersJaccard<-function(d) { 
  
  tf1<-chea1[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  
  tf2<-chea1[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  
  jaccard<-(1-(length(intersect(tf1,tf2))/length(union(tf1,tf2))))
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  row<-data.frame(d[1],d[2],pVal[[1]],jaccard,length(genes),length(tf1),length(tf2),length(intersect(tf1,tf2)),stringsAsFactors = FALSE)
  colnames(row)=c("tf1","tf2","pVal","jaccard","background genes","tf1 target genes","tf2 target genes","intersection size")
  return(row)
}

cheaPVal<-as.data.frame(do.call(rbind,(lapply(cheaColsCombo, function(d) fishersJaccard(d))))) # PRACTICE

proc.time() - ptm

write.table(cheaPVal,"~/cheaPVal.tsv",col.names=TRUE,sep="\t")

cheaPVal$match<-NA
cheaPVal[["match"]]<-apply(cheaPVal,1,function(x) {ifelse(gsub("_.*","",x[1])==gsub("_.*","",x[2]),1,0)})
write.table(cheaPVal,"~/cheaPValJMatch.tsv",col.names=TRUE,sep="\t")

## Everything that was in fishersExact.R before I renamed creedsFisher


dataset1=chea1
dataset2=encode

cheaColsCombo<-combn(colnames(chea1),2,simplify=FALSE)
encodeColsCombo<-combn(colnames(encode),2,simplify=FALSE)

## chea and creeds
bothColsCombo<-combn(union(colnames(dataset1),colnames(dataset2)),2,simplify=FALSE)
bothColsCombo<-setdiff(bothColsCombo,dataset1)
bothColsCombo<-setdiff(bothColsCombo,dataset2)

genes<-unique(union(unlist(dataset1),unlist(dataset2)))%>%na.omit()
c<-bothColsCombo[1:100]

b<-lapply(bothColsCombo,function(x){
  
  if(any(grepl(x[1],colnames(dataset2)))){
    #x<-rev(x)
    return(x)
  }
  #return(x)
})


## encode and creeds
# bothColsCombo<-combn(union(colnames(encode),colnames(creedsFisher)),2,simplify=FALSE)
# bothColsCombo<-setdiff(bothColsCombo,encodeColsCombo)
# bothColsCombo<-setdiff(bothColsCombo,creedsColsCombo)
# 
# genes<-unique(union(unlist(encode),unlist(creedsFisher)))%>%na.omit()

## chea and encode
bothColsCombo<-combn(union(colnames(chea1),colnames(encode)),2,simplify=FALSE)
bothColsCombo<-setdiff(bothColsCombo,cheaColsCombo)
bothColsCombo<-setdiff(bothColsCombo,encodeColsCombo)

genes<-unique(union(unlist(chea1),unlist(encode)))%>%na.omit()

c<-bothColsCombo[1:100]
##
ptm <- proc.time()

b<-lapply(bothColsCombo,function(x){
  
  if(any(grepl(x[1],colnames(encode)))){
    #x<-rev(x)
    return(x)
  }
  #return(x)
})
proc.time() - ptm
##
#bothPVal<-data.frame(tf1=character(),tf2=character(),pVal=character(),jaccard=character(),stringsAsFactors = FALSE)
ptm <- proc.time()

fishersJaccard<-function(d) { # will have to change function for comparing chea and encode
  
  tf1<-chea1[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  
  tf2<-encode[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  
  jaccard<-(1-(length(intersect(tf1,tf2))/length(union(tf1,tf2))))
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  row<-data.frame(d[1],d[2],pVal[[1]],jaccard,stringsAsFactors = FALSE)
  colnames(row)=c("tf1","tf2","pVal","jaccard")
  return(row)
}

bothPVal<-as.data.frame(do.call(rbind,(lapply(bothColsCombo, function(d) fishersJaccard(d)))))

proc.time() - ptm

write.table(bothPVal,"~/cheaEncodePVal.tsv",col.names=TRUE,sep="\t")
#encodePVal<-read.table("~/encodePVal.tsv",header=TRUE,stringsAsFactors = FALSE)

## Making the match column
bothPVal$match<-NA
bothPVal[["match"]]<-apply(bothPVal,1,function(x) {ifelse(gsub("_.*","",x[1])==gsub("_.*","",x[2]),1,0)})
# bothPVal[["match"]]<-apply(bothPVal,1,function(x) {ifelse(gsub("_.*","",x[1])==gsub(".* ","",x[2]),1,0)})

write.table(bothPVal,"~/cheaEncodePValJMatch.tsv",col.names=TRUE,sep="\t")

cheaCreedsPVal<-bothPVal
encodeCreedsPVal<-bothPVal
cheaEncodePVal<-bothPVal

## 
a<-apply(cheaCoexp,2,rocAUC,term_list=cheaGenes,score_list_labels=rownames(cheaCoexp),top_rank="greatest")

b<-apply(encodeCoexp,2,rocAUC,term_list=encodeGenes,score_list_labels=rownames(encodeCoexp),top_rank="greatest")

par(mfrow = c(2,1))
a_FET = cheaPVal[order(cheaPVal$pVal),]
a_FET$rankForPlot = 1:nrow(a_FET)
b<-ggplot(a_FET, aes(x=rankForPlot, group = as.character(match), color = as.character(match))) + geom_density() + xlab("plot1")
dataset_FET = encodePVal[order(encodePVal$pVal),]
dataset_FET$rankForPlot = 1:nrow(dataset_FET)
c<-ggplot(dataset_FET, aes(x=rankForPlot, group = as.character(match), color = as.character(match))) + geom_density() + xlab("plot2")

sample<-list(b,c)

#plotsList<-lapply(datasetsList,plots)
#grid.arrange(plotsList[[1]],plotsList[[2]])

# ggplot(creedsEncode_FET, aes(x=rank, group = as.character(match), color = as.character(match))) + geom_density()
# 
# ggplot(dataset_FET, aes(x=rank/332520, group = as.character(match), color = as.character(match))) + geom_density()

# encodeCoexpAUC<-read.table("~/encodeCoexpAUCMatch.tsv",header=T,stringsAsFactors = F)

#+ geom_text(aes(label=text),family="Helvetica")
#+ theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9),axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),legend.text = element_text(size=7))  

# 
# pValSimple<-function(tf1,tf2) {
#   print(head(tf1))
#   contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
#   return(fisher.test(contTable,alternative = "two.sided",conf.int=FALSE))
# }
# 
# a<-apply(creedsFisher[,1:2],2,function(tf1) apply(datasetIntersection[,1:2],2,function(tf2)pValSimple(tf2),tf1=tf1))
#   
# a<-apply(datasetIntersection,2,function(tf2)pValSimple,tf1=creedsFisher[,1])
# 
# a<-apply(datasetIntersection,2,function(tf2){return(length(intersect(tf1,tf2)))},tf1=creedsFisher[,1])
# 
# a<-apply(datasetIntersection,2,function(tf2){
#   apple<-intersect(tf2,tf1)
#   return(length(apple))},tf1=c("a","b"))

## Code from https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters ; apply function with multiple parameters
# mylist <- list(a=1,b=2,c=3)
# myfxn <- function(var1,var2){
#   var1*var2
# }
# 
# a<-sapply(mylist,myfxn,var2=2)

## read in pVal files
cheaTest<-read.table("~/cheaTest.tsv",header=T)
encodeTest<-read.table("~/encodeTest.tsv",header=T)
coexpTest<-read.table("~/coexpTest.tsv",header=T)

## intersection
intersection<-intersect(colnames(coexp),colnames(chea2)) ## ignoring TCF12.HEB
intersection<-intersect(one,colnames(encode2))

intersection<-gsub("\\..*","",intersection)
intersection<-(unique(intersection))

# rownames(cheaTest1)=gsub("\\..*","",rownames(cheaTest1))
c<-colnames(chea2)[colnames(chea2)%in%colnames(coexp)]
e<-unique(c)
d<-colnames(chea2)%in%colnames(coexp)

f<-colnames(chea2)[intersection==T]
g<-unique(f)

g<-sort(g)
e<-sort(e)
length(e)=length(g)
apple<-cbind(e,g)
apple<-as.data.frame(apple,stringsAsFactors = F)

h<-colnames(coexp)
h<-sort(h)
length(g)=length(h)
banana<-cbind(g,h)

datasetTest<-cheaTest1

j<-gsub("\\..*","",colnames(datasetIntersection))

## This code gets the missing ones (TCF12.HEB) in ChEA
dataset=chea1
intersection<-sapply(colnames(dataset),function(x){return(gsub("\\..*","",x)%in%colnames(coexp))}) # would have to sapply to colnames(gsub(other dataset)) for real intersection
datasetIntersection<-dataset[,intersection==T]

##
pValSimple(creedsFisher[,2],datasetIntersection[,2])
# b<-sapply(intersection[1],sample)
# sample<-function(x){
#   return(min(x,x+1))
# }



## later: cheaTest1 (save as another one)
## rownames(datasetTest)=as.vector(gsub("\\..*","",rownames(datasetTest)))
## duplicated[datasetTest] to find the duplicated
# sapply(duplicated,function(x) return(max(datasetTest[grep(x,rownames(datasetTest),value=T),])))

rankedSets<-lapply(rankedSets,subsets) ## doesn't actually change the things inside the list. Just changes this ranked set collection.

subsets<-function(x) {
  return(x[,colnames(x)%in%colnames(creedsIntersection)])
}

## add columns (creedsFisher colnames) to datasets
datasetsList2<-lapply(datasetsList2,function(x){
  colnames(x)=colnames(creedsFisher)
  return(x)
})

#  
# print(rownames(as.data.frame(dataset1))[a])
# print(rownames(as.data.frame(dataset2))[1:3])
#print(b)
# print(cat("c is ",c))
# if(identical(c,integer(0))) {
#   d<-(as.data.frame(newDataset))[a,b]
#   #print(cat("d is ",d))
#   return(d)
# }
# else {
#   e<-min((as.data.frame(dataset1))[a,b],(as.data.frame(dataset2))[c,b])
#   #print(cat("e is",e))
#   return(e)
# }

# for(i in 1:length(datasetsList3)) {
#   # minPVal1<-as.data.frame(list[[1]])
#   # minPVal2<-as.data.frame(list[[2]])
#   # toAdd<-minPVal(minPVal1,minPVal2)
#   print(dim(datasetsList[[1]]))
# }

# x<-minPVal(g,h)

g<-matrix(c(0,0,0,0),nrow=2)
h<-matrix(c(1,1,1,1,1,1),nrow=3)
i<-matrix(c(2,2,2,2,2,2),nrow=3)
j<-list(g,h,i)

k<-combn(j,2,sample)

sample<-function(x) {
  print(x[[1]])
  # print("next")
  # print(typeof(x[1]))
  print(x[[2]])
  # print(typeof(x[2]))
  print("break")
}

e<-as.vector(unlist(rowNames[1]))

d<-lapply(e,function(a) {
  # print(a)
  # print("next")
  return(a)})

## get one column based on another? from stackoverflow
result <- data.frame(a, b=(ifelse(a=="foo","found",b)))

#datasetsPlot1<-list(Combined_Datasets=ranksCombinedDatasetPlot,Greatest_PValue=ranksGreatestPValPlot,Average_PValue=ranksAveragePValPlot,Mean_Ranks_Not_Normalized=ranksMeanNotNormalizedPlot,Mean_Ranks_Normalized=ranksMeanNormalizedPlot,Lowest_PValue=ranksLowestPValPlot,Smallest_Dataset=ranksSmallestDatasetPlot,Multiplied_PValue=ranksMultipliedPValPlot,ranksLowestRankPlot)

p1<-ggplot(data=datasetsPlot1[[1]],aes(x = value, color = variable)) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Combined Datasets",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5))+ theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold",windowsFont("Verdana"))) + geom_line(aes(color=variable), stat="density", size=0.705, alpha=0.8) #,expand = c(0,0) + geom_density(size=0.705,alpha=0.1) 
p2<-ggplot(data=datasetsPlot1[[2]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Greatest P-Value",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p3<-ggplot(data=datasetsPlot1[[3]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Average P-Value",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5))  + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold"))#,expand = c(0,0)
p4<-ggplot(data=datasetsPlot1[[4]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Mean Ranks (Not Scaled)",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p5<-ggplot(data=datasetsPlot1[[5]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Mean Ranks (Scaled)",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p6<-ggplot(data=datasetsPlot1[[6]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Lowest P-Value",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p7<-ggplot(data=datasetsPlot1[[7]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Smallest Dataset",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p8<-ggplot(data=datasetsPlot1[[8]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Multiplied P-Values",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)
p9<-ggplot(data=datasetsPlot1[[9]],aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Lowest Rank",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold")) #,expand = c(0,0)

##

t<-t(sapply(1:length(colnames(as.data.frame(ranks))),example))

example<-function(x){
  r<-ranks[x,x]
  s<-2*r
  return(data.frame(exp=r,control=s))
}

## old normalize function
normalize<-function(x) {
  largest<-max(x)
  return(sapply(x,function(y) return(y/largest)))
}

##
# colnames(datasetsList1[[1]])<-colnames(creedsFisher) despite it not saying so I think the column names really do have dashes.
# colnames(f)<-colnames(creedsFisher)

# write.table(f,"~/datasetsList1-Coexp.tsv",col.names=T)

# datasetsList1<-lapply(datasetsList1,function(x){
#   colnames(x)=colnames(creedsFisher)
#   return(x)
# })

##
## normalize Test whether need or not
datasetsList4<-lapply(datasetsList3,function(a)apply(a,2,normalize))

normalize<-function(x) {
  largest<-max(x)
  return(sapply(x,function(y) return(y/largest)))
}

j<-as.data.frame(datasetsList4[[1]])
k<-as.data.frame(datasetsList4[[2]])
l<-as.data.frame(datasetsList4[[3]])

## 
## compare
m<-as.data.frame(sample[[1]])
n<-as.data.frame(sample[[2]])
o<-as.data.frame(sample[[3]])

datasetsList5<-list()
sample<-combn(datasetsList3,2,wrapper)
datasetsList5<-sample

# write.table(three,"~/datasetsListMeanRanks-Three.tsv",col.names=T)

wrapper<-function(list) {
  minPVal1<-as.data.frame(list[[1]])
  minPVal2<-as.data.frame(list[[2]])
  toAdd<-minRank(minPVal1,minPVal2)
  datasetsList5[[length(datasetsList5)+1]]<-toAdd
  return(datasetsList5)
}

minRank<-function(dataset1,dataset2) {
  
  rows1<-gsub("-.*","",rownames(as.data.frame(dataset1)))
  rows2<-gsub("-.*","",rownames(as.data.frame(dataset2)))
  together<-as.vector(sort(union(rows1,rows2)))
  
  newDataset<-data.frame(matrix(ncol=length(colnames(as.data.frame(dataset1))),nrow=length(together)),row.names=together)
  
  test<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset1))),function(d) sapply(1:length(rownames(as.data.frame(newDataset))),minBetweenDatasets,b=d,dataset1=dataset1,dataset2=dataset2,newDataset=newDataset))))) #length(colnames(as.data.frame(dataset1))) #length(rownames(as.data.frame(newDataset)))
  colnames(test)<-colnames(dataset1)
  rownames(test)<-rownames(newDataset)
  return(test)
}

minBetweenDatasets <-function(a,b,dataset1,dataset2,newDataset) {
  c<-paste0("^",rownames(as.data.frame(newDataset))[a],"$")
  first<-grep(c,rownames(as.data.frame(dataset1)))
  second<-grep(c,rownames(as.data.frame(dataset2)))
  
  if(identical(first,integer(0))&&!identical(second,integer(0))) {
    return((as.data.frame(dataset2))[second,b])
  } else if(identical(second,integer(0))&&!identical(first,integer(0))) {
    return((as.data.frame(dataset1))[first,b])
  } else {
    return(min((as.data.frame(dataset1))[first,b],(as.data.frame(dataset2))[second,b])) #((as.data.frame(dataset1))[first,b]+(as.data.frame(dataset2))[second,b])/2 #((as.data.frame(dataset1))[first,b]+(as.data.frame(dataset2))[second,b])/2
  }
  
}

three<-minRank(as.data.frame(datasetsList5[[1]]),as.data.frame(datasetsList3[[3]])) # have to do something different for three: (as.data.frame(dataset1))[first,b]+(((2*(as.data.frame(dataset2))[second,b])-(2*(as.data.frame(dataset1))[first,b]))/6) If change normalized ranks to not, have to change from datasetsList4 to 3
datasetsList5[[length(datasetsList5)+1]]<-three
datasetsList6<-list(as.data.frame(datasetsList3[[1]]),as.data.frame(datasetsList3[[2]]),as.data.frame(datasetsList3[[3]]),as.data.frame(datasetsList5[[1]]),as.data.frame(datasetsList5[[2]]),as.data.frame(datasetsList5[[3]]),as.data.frame(datasetsList5[[4]])) # may have to change datasetsList4 to 3

### datasetsList3<-lapply(datasetsList2,function(x) apply(x,c(1,2),as.integer))

#
if(identical(first,integer(0))&&(!identical(second,integer(0))&&!identical(third,integer(0)))) {
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

##
minBetweenDatasets<-function(b,together,dataset1,dataset2,maximum) {
  c<-paste0("^",together[b],"$")
  
  first<-grep(c,colnames(as.data.frame(dataset1)))
  second<-grep(c,colnames(as.data.frame(dataset2)))
  
  combined<-rep(NA,maximum)
  
  # print(b)
  
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
      return(combined)
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
      return(combined)
    }
  }
  if(!identical(first,integer(0))&&!identical(second,integer(0))) {
    # print(c)
    # print(piece1[1:5]) # fix, throw it out
    # print(piece2[1:5])
    
    piece<-unique(intersect(piece1,piece2))
    # print(piece)
    # print("-----")
    if(identical(piece,character(0))) {
      return(combined)
    }
    else {
      combined[1:length(piece)]<-piece
      # print(length(combined))
      return(combined)
    }
  }
}
##
minBetweenDatasets<-function(b,together,dataset1,dataset2,maximum) {
  c<-paste0("^",together[b],"$")
  
  first<-grep(c,colnames(as.data.frame(dataset1)))
  second<-grep(c,colnames(as.data.frame(dataset2)))
  
  combined<-rep(NA,maximum)
  
  # print(b)
  
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
      return(combined)
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
      return(combined)
    }
  }
  if(!identical(first,integer(0))&&!identical(second,integer(0))) {
    # print(c)
    # print(piece1[1:5]) # fix, throw it out
    # print(piece2[1:5])
    
    piece<-unique(intersect(piece1,piece2))
    # print(piece)
    # print("-----")
    if(identical(piece,character(0))) {
      return(combined)
    }
    else {
      combined[1:length(piece)]<-piece
      # print(length(combined))
      return(combined)
    }
  }
}

##
minBetweenDatasets<-function(b,together,dataset1,dataset2,maximum) {
  c<-paste0("^",together[b],"$")
  
  first<-grep(c,colnames(as.data.frame(dataset1)))
  second<-grep(c,colnames(as.data.frame(dataset2)))
  
  combined<-rep(NA,maximum)
  
  # print(b)
  
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
    # print(c)
    # print(piece1[1:5]) # fix, throw it out
    # print(piece2[1:5])
    
    piece<-unique(intersect(piece1,piece2))
    # print(piece)
    # print("-----")
    if(identical(piece,character(0))) {
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
    else {
      combined[1:length(piece)]<-piece
      # print(length(combined))
      temp[[length(temp)+1]]<-combined
      return(temp)
    }
  }
}

plots2<-function(dataset,name) {
  
  ggplot(data=dataset,aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title=name,x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5)) + theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold"))+coord_cartesian(expand = c(0,0))
}

datasetsPlot1<-list("Top Rank Scaled"=ranksTopRankPlotScaled,"Top Rank Unscaled"=ranksTopRankPlotUnscaled)
plotsList1<-lapply(seq_along(datasetsPlot1),function(i){plots1(datasetsPlot1[[i]],names(datasetsPlot1)[i])})
p1<-ggplot(data=subset(ranksTopRankPlotUnscaled,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + 
  scale_colour_pander() + labs(color="Method") + labs(title="Top Rank Unscaled",x="Rank of TF",y="Density") + 
  coord_cartesian(xlim = c(0, 420), ylim = c(0,0.01),expand = c(0,0)) + 
  theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),
        plot.title = element_text(size=11,hjust = 0,face="bold"))

p2<-ggplot(data=subset(ranksTopRankPlotScaled,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + 
  scale_colour_pander() + labs(color="Datasets") + labs(title="Top Rank Scaled",x="Rank of TF",y="Density") + 
  coord_cartesian(xlim = c(0, 1.1), ylim = c(0,4)) + 
  theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),
        plot.title = element_text(size=11,hjust = 0,face="bold"))+coord_cartesian(expand = c(0,0)) 

grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)

p1<-ggplot(data=subset(ranksTopRankPlotUnscaled,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + scale_colour_pander() + labs(color="Method") + labs(title="Mean Rank Unscaled",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 420), ylim = c(0,0.01),expand = c(0,0)) + theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold"))
p2<-ggplot(data=subset(ranksMeanRankPlotScaled,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.705) + theme_classic() + scale_colour_pander() + labs(color="Datasets") + labs(title="Mean Rank Scaled",x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,4)) + theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),plot.title = element_text(size=11,hjust = 0,face="bold"))+coord_cartesian(expand = c(0,0)) 
grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)

## other normalize function
normalize<-function(x,maximum) {
  x<-as.numeric(x)
  return(x/maximum)
}

#
## rank again
datasetsList5<-lapply(datasetsList5,function(x) apply(x,2,rank))

j<-as.data.frame(datasetsList5[[1]])
k<-as.data.frame(datasetsList5[[2]])
l<-as.data.frame(datasetsList5[[7]])

## combining them
ranks<-read.table("~/ranksLowestRank.tsv",header=T)

ranks<-as.data.frame(do.call(cbind,(lapply(datasetsList5, function(t) t(sapply(1:length(colnames(as.data.frame(t))),ranked,d=t)))))) #length(colnames(as.data.frame(t)))
rownames(ranks)=colnames(creedsFisher)
colnames(ranks)=c("ChEA","ChEA Control","ENCODE","ENCODE Control","Coexp","Coexp Control","ChEA Encode","ChEA ENCODE Control","ChEA Coexp","ChEA Coexp Control","ENCODE Coexp","ENCODE Coexp Control","All three","All three Control")
ranks<-apply(ranks,2,as.numeric)
ranks<-as.data.frame(ranks,stringsAsFactors = F)
ranks<-as.matrix(ranks)

write.table(ranks,"~/ranksIntersection.tsv",col.names = T,row.names=T)

ranked<-function(x,d){ 
  a<-which(rownames(as.data.frame(d))==gsub("-.*","",colnames(as.data.frame(d))[x]))
  # print(paste("a is",a))
  if(identical(a,integer(0))) {
    # print(max(d[,x]))
    r<-NA
  }
  else {
    # print(d[a,x])
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  # print(paste(r,s))
  return(data.frame(exp=r,control=s))
}

ranksLowestPVal<-read.table("~/ranksLowestPVal.tsv",header=T,stringsAsFactors=F)
ranksCombinedDataset<-read.table("~/ranksCombinedDatasets.tsv",header=T,stringsAsFactors=F)

## Visualize data
library(reshape2)

# ranks=ranks[,c("ChEA","ENCODE","Coexp","ChEA ENCODE","ChEA Coexp","ENCODE Coexp","All three","All three Control","ChEA Control","ENCODE Control","Coexp Control","ChEA ENCODE Control","ChEA Coexp Control","ENCODE Coexp Control")]
ranks=ranks[,c(1,3,5,7,9,11,13,14,2,4,6,8,10,12)]

maximum<-max(as.vector(unlist(ranks)),na.rm = T)
ranksNormalize<-apply(ranks,c(1,2),normalize,maximum=maximum)
ranksNormalize<-as.data.frame(ranksNormalize,stringsAsFactors = F)

normalize<-function(x,maximum) {
  x<-as.numeric(x)
  return(x/maximum)
}

ranksMeltNormalize<-melt(ranksNormalize)

library(ggplot2)
library(ggthemes)

plot_subset = c("ChEA","ChEA Control", "ChEA Coexp", "Coexp Control", "Coexp")
subset_ranksMeltScaled = ranksMeltNormalize[is.element(ranksMeltNormalize$variable,plot_subset),]

ggplot(data=subset(ranksMeltNormalize,!is.na(value)),aes(x = value, color = variable)) + geom_density(size=0.75) + theme_classic() + scale_colour_pander() + labs(color="Method") + labs(x="Rank of TF",y="Density") + coord_cartesian(xlim = c(0, 1.1), ylim = c(0,3.5),expand = c(0,0)) 

hist(ranksLargeNormalize$ENCODE)

##
# make_union <-function(b,together,dataset1,dataset2,maximum) {
#   c<-paste0("^",together[b],"$")
# 
#   print(together[b])
# 
#   first<-grep(c,colnames(as.data.frame(dataset1)))
#   second<-grep(c,colnames(as.data.frame(dataset2)))
# 
#   combined<-rep(NA,maximum)
#   piece=character()
#   if(!identical(first,integer(0))) {
#     # print("1")
#     # extracted1<-character()
#     # extracted1<-as.data.frame(do.call(rbind,(sapply(first,extract,dataset=as.data.frame(dataset1)))))
#     piece<-sapply(first,extract,dataset=as.data.frame(dataset1),piece)
#     # print(extracted1[1:5])
#     # print(extracted1[1:5])
#     # extracted1<-as.data.frame(do.call(rbind,(sapply(first,extract,dataset=as.data.frame(dataset1)))))
#     # print(paste("extract dim",dim(extracted1)))
#     # combined[1:length(extracted1)]<-extracted1[1:5,1]
# 
#   }
#   if(!identical(second,integer(0))) {
#     # print("2")
#     # extracted2<-as.data.frame(do.call(rbind,(sapply(second,extract,dataset=as.data.frame(dataset2)))))
#     piece<-sapply(second,extract,dataset=as.data.frame(dataset2),piece)
#     # if(identical(first,integer(0))) {
#     #   combined[1:length(rownames(extracted2))]<-extracted2[,1]
#     # }
#     # else{
#     #   combined[(length(rownames(extracted1))+1):(length(rownames(extracted1))+length(rownames(extracted2)))]<-extracted2[,1]
#     # }
# 
#   }
#   combined[1:length(piece)]<-piece
#   return(combined)
# }
# 
# extract<-function(x,dataset,piece) {
#   noNA<-as.data.frame(dataset)[,x]
#   noNA<-noNA[!is.na(noNA)]
#   noNA<-as.vector(noNA)
#   # p<-data.frame("tf"=noNA)
#   # vector<-rbind(vector,noNA)
#   piece=c(piece,noNA)
#   # print(paste("p",p[1:5,1]))
#   # print(dim(p))
#   return(piece)
# }

##
write.table(chea2,"~/cheaEdited.tsv",sep="\t",col.names=TRUE)

##
# andNetworkCols<-colnames(andNetwork)
# write.table(andNetworkCols,"C:/Users/maayanlab1/Downloads/andNetworkCols.tsv",sep="\t",row.names = FALSE,col.names = FALSE)
# andNetworkRows<-rownames(andNetwork)
# write.table(andNetworkRows,"C:/Users/maayanlab1/Downloads/andNetworkRows.tsv",sep="\t",row.names = FALSE,col.names = FALSE)
# c<-read.table("C:/Users/maayanlab1/Downloads/andNetwork.tsv",sep="\t", row.names=andNetworkRows,col.names = andNetworkCols, stringsAsFactors = FALSE)
# c<-read.csv("C:/Users/maayanlab1/Downloads/andNetwork.csv",row.names=1)

#pearsonCorr<-cor(expression)
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


##
## get supplementary data
write(tfs,"~/tfs.txt",sep="\t")
write(colnames(chea2),"~/colnames(chea2).txt",sep="\t")
write(cheaGenes,"~/cheaGenes.txt",sep="\t")
write(colnames(encode2),"~/colnames(encode2).txt",sep="\t")
write(encodeGenes,"~/encodeGenes.txt",sep="\t")

##


sapply(colnames(chea1)[1:5], function(x) {ifelse(gsub("_.*","",colnames(chea1)[1])==gsub("_.*","",x),print(colnames(chea[x])),print(""))})

length(colnames(chea1))
length(colnames(encode))
a<-union(colnames(chea1),colnames(encode))

##ranks_subset<- ranksSmallIntersectionPlotScaled %>% filter(variable %in% c("ENCODE Coexp","ChEA Coexp")==T)

### insertList<-list()
# insertList[[1]]<-datasetsListCombine2[[1]]
# insertList[[2]]<-datasetsListCombine[[3]]

##
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

##
for(i in 0:4) {
  print(sum(!is.na(ranks[,2*i+1])))
}

exp_distribution_df_1<-exp_distribution_df #[rev(order(exp_distribution_df$total)),]

for(i in 2:length(datasetsList2_ppi)) {
  a<-combn(JUN,i,function(x)Reduce(intersect,x))
  jun[[length(jun)+1]]<-a
}

jun_df<-data.frame(sets=character(),stringsAsFactors = F) #not intersection yet

jun_labels<-list("chea1","chea2","encode1","encode2","encode3","encode4","encode5","encode6","coexp1","bioplex1","humap1","oldcreeds1","oldcreeds2","oldcreeds3","oldcreeds4")
stat3_labels<-list("chea1","chea2","chea2","chea2","chea2","encode4","encode5","encode6","coexp1","bioplex1","humap1","oldcreeds1","oldcreeds2","oldcreeds3","oldcreeds4")

tf_exps_intersect <- function(list) {
  unlisted<-unlist(list)
  jun_df = rbind(jun_df,paste(unlisted[1],unlisted[2],sep="-"))
  return(jun_df)
}

jun_df<-combn(jun_labels,2,tf_exps_intersect)

jun_df1<-jun_df
jun_df1<-unlist(jun_df)
jun_df1<-as.vector(jun_df1)
jun_df1<-as.data.frame(jun_df1,stringsAsFactors = F)

# method_df<-function(dataset1,dataset2,method) {
#   
#   rows1<-gsub("-.*","",rownames(as.data.frame(dataset1)))
#   rows2<-gsub("-.*","",rownames(as.data.frame(dataset2)))
#   together<-as.vector(sort(union(rows1,rows2)))
#   
#   test<-as.data.frame(do.call(cbind,(lapply(1:length(colnames(as.data.frame(dataset1))),function(d) sapply(1:length(together),method_compute,b=d,dataset1=dataset1,dataset2=dataset2,together=together,method=method))))) #length(colnames(as.data.frame(dataset1))) #length(rownames(as.data.frame(newDataset)))
#   colnames(test)<-colnames(dataset1)
#   rownames(test)<-together
#   return(test)
# }

ranked<-function(x,d){ 
  b<-gsub("-.*","",colnames(as.data.frame(d))[x])
  a<-grep(paste0("^",b,"$"),rownames(as.data.frame(d)),ignore.case = T)
  # a<-which(rownames(as.data.frame(d))==gsub("-.*","",colnames(as.data.frame(d))[x]))
  # print(b)
  # print(a)
  # print(paste0("x is",x))
  # print(rownames(as.data.frame(d))[1:5])
  if(identical(a,integer(0))) {
    r<-NA
  }
  else {
    r<-d[a,x]
  }
  s<-sample(as.data.frame(d)[,x],1)
  return(data.frame(exp=r,control=s))
}

## Trying to do it individually
bioplex_pvals_old_creeds<-makepVals(bioplex,creedsFisher) # because you wrote over the other bioplex_pvals

ptm<-proc.time()
try_2<-makepVals(humap,new_creeds)
proc.time()-ptm

write.table(new_creeds,"~/new_creeds.tsv",sep="\t",col.names = T,row.names = T)

try_2<-read.table("~/humap_ppi_new_creeds.tsv",sep="\t",header=T,stringsAsFactors = F)

# c<-toupper(gsub("\\.","-",colnames(humap_pval)))
# colnames(humap_pval)<-c
# colnames(bioplex_pval)<-c
# colnames(datasetsList0_ppi_new_creeds[[3]])<-c

datasetsListMultipliedRanks_ppi_2<-list(datasetsListMultipliedPVal[[1]],datasetsListMultipliedPVal[[2]],datasetsListMultipliedPVal[[3]],
                                        datasetsListMultipliedPVal[[4]],datasetsListMultipliedPVal[[5]],datasetsListMultipliedPVal[[6]],
                                        datasetsListMultipliedPVal[[7]],datasetsListMultipliedPVal[[8]],datasetsListMultipliedPVal[[9]],
                                        datasetsListMultipliedPVal[[10]])

datasetsListMultipliedRanks_ppi_2_ranked<-list()

for(i in 1:10) {
  datasetsListMultipliedRanks_ppi_2_ranked[[i]]<-apply(datasetsListMultipliedRanks_ppi_2[[i]],2,rank)
}

a<-gsub("\\.","-",colnames(as.data.frame(datasetsList6_ppi_new_creeds[[5]])))
colnames(datasetsList6_ppi_new_creeds[[5]])<-a

ptm<-proc.time()
sample<-t(sapply(1:length(colnames(e)),ranked,d=e))
proc.time()-ptm

sample<-as.matrix(sample)

sample<-t(sapply(1:length(colnames(e)),ranked,d=e))
rownames(sample)<-a

humap_rank<-sample
humap_rank<-as.data.frame(humap_rank)

bioplex_rank<-sample
bioplex_rank<-as.matrix(bioplex_rank)
bioplex_rank<-as.data.frame(bioplex_rank)

coexp_rank<-sample
coexp_rank<-as.matrix(coexp_rank)
coexp_rank<-as.data.frame(coexp_rank)

encode_rank<-sample
encode_rank<-as.matrix(encode_rank)
encode_rank<-as.data.frame(encode_rank)

chea_rank<-sample
chea_rank<-as.matrix(chea_rank)
chea_rank<-as.data.frame(chea_rank)



