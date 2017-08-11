install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggthemes")
install.packages("extrafont")
install.packages("Grid")

library(ggplot2)
library(gridExtra)
library(ggthemes)
library(extrafont)
library(extrafontdb)
library(grid)

# font_import()
#loadfonts(device="postscript")

plots<-function(dataset,name) {
  
  dataset_FET = dataset[order(dataset$pVal),]
  dataset_FET$rankForPlot = 1:nrow(dataset_FET)
  dataset_FET$rankForPlot<-normalize(dataset_FET$rankForPlot)
  
  ggplot(dataset_FET, aes(x=rankForPlot, group = as.character(match), color = as.character(match))) + geom_density(size=0.73) + 
    labs(title=name,x="P-Value Rank (Fisher's Exact Test)",y="Density") + theme_classic() + scale_colour_pander() + labs(color="Match") + 
    theme(axis.text.x = element_text(size=8,hjust=1),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),
          axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),
          plot.title = element_text(size=11,hjust = 0,face="bold"))+coord_cartesian(expand = c(0,0)) 
  #xlim = c(0, 1.1), ylim = c(0,.00005), #theme(axis.text.x = element_text(size=8,angle=90,hjust=1)    
}

datasetsPlot<-list("ENCODE Co-expression"=encodeCoexp,"ENCODE"=encodePVal,"ChEA Co-expression"=cheaCoexp,"ChEA"=cheaPVal,
                   "CREEDS Co-expression"=creedsCoexp,"CREEDS"=creedsPVal,"ChEA ENCODE"=cheaEncodePVal,"ChEA CREEDS"=cheaCreedsPVal,
                   "ENCODE CREEDS"=encodeCreedsPVal)

plotsList<-lapply(seq_along(datasetsPlot),function(i){plots(datasetsPlot[[i]],names(datasetsPlot)[i])})
# grid.arrange(grobs=plotsList)
grid_arrange_shared_legend(plotsList[[1]], plotsList[[2]], plotsList[[3]], plotsList[[4]], plotsList[[5]], plotsList[[6]], plotsList[[7]], plotsList[[8]], plotsList[[9]], ncol = 3, nrow = 3)

#---------
plots1<-function(dataset,name) {
  
  ggplot(data=subset(dataset,!is.na(value)),aes(x = value, color = variable, linetype = variable)) + geom_density(size=0.705) + theme_classic() + 
    labs(color="Datasets",linetype = "Datasets") + labs(title=name,x="Rank of TF",y="Density") + 
    coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5),expand = c(0,0)) + 
    theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),
      axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),
      plot.title = element_text(size=11,hjust = 0,face="bold")) +
    scale_linetype_manual(values=c("solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted")) +
    scale_colour_manual(values = c("deepskyblue","deepskyblue","seagreen","seagreen","yellow","yellow","navy","navy","red3","red3","hotpink","hotpink","grey58","grey58"))
}

datasetsPlot1<-list("Target Genes Intersection"=ranksIntersectionPlotScaled,"Mean Ranks"=ranksMeanRankPlotScaled,
                    "Top Rank"=ranksTopRankPlotScaled,"Multiplied P-Values"=ranksMultipliedPValPlotScaled)
plotsList1<-lapply(seq_along(datasetsPlot1),function(i){plots1(datasetsPlot1[[i]],names(datasetsPlot1)[i])})
grid_arrange_shared_legend(plotsList1[[1]], plotsList1[[2]], plotsList1[[3]], plotsList1[[4]], ncol = 2, nrow = 2)

#---------
plots2<-function(dataset,name) {
  
  ggplot(data=subset(dataset,!is.na(value)),aes(x = value, color = variable, linetype = variable)) + geom_density(size=0.705) + theme_classic() + 
    labs(color="Datasets",linetype = "Datasets") + labs(title=name,x="Rank of TF",y="Density") + 
    coord_cartesian(xlim = c(0, 1.1), ylim = c(0,2.5),expand = c(0,0)) + 
    theme(axis.text.x = element_text(size=8,hjust=0),axis.text.y = element_text(size=8),axis.title.x = element_text(size=9),
          axis.title.y = element_text(size=9),legend.title = element_text(size=9),legend.text = element_text(size=8),
          plot.title = element_text(size=11,hjust = 0,face="bold")) +
    scale_linetype_manual(values=c("solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted")) +
    scale_colour_manual(values = c("deepskyblue","deepskyblue","seagreen","seagreen","yellow","yellow","navy","navy","red3","red3","hotpink","hotpink","grey58","grey58"))
}

datasetsPlot2<-list("Target Genes Intersection"=ranksSmallIntersectionPlotScaled,"Mean Ranks"=ranksSmallMeanRankPlotScaled,"Top Rank"=ranksSmallTopRankPlotScaled,"Multiplied P-Values"=ranksSmallMultPValPlotScaled)
plotsList2<-lapply(seq_along(datasetsPlot2),function(i){plots2(datasetsPlot2[[i]],names(datasetsPlot2)[i])})
grid_arrange_shared_legend(plotsList2[[1]], plotsList2[[2]], plotsList2[[3]], plotsList2[[4]], ncol = 2, nrow = 2)

#--------- for changing the datasets
ranks<-read.table("~/ranksLarge.tsv",header=T)

ranks$Random.Rank = sample(length(rownames(ranks)))
colnames(ranks)=c("ChEA","ENCODE","Coexp","ChEA Encode","ChEA Coexp","ENCODE Coexp","All three","Random Ranks")

maximum<-max(as.vector(unlist(ranks)),na.rm = T)
ranksNormalize<-apply(ranks,c(1,2),normalize,maximum)
ranksNormalize<-as.data.frame(ranksNormalize,stringsAsFactors = F)

normalize<-function(x,maximum) {
  x<-as.numeric(x)
  return(x/maximum)
}

ranksMeltNormalize<-melt(ranksNormalize)

ranksMultipliedPValPlot<-ranksMeltNormalize

write.table(ranksMultipliedPValPlot,"~/ranksMultipliedPValPlot.tsv",col.names = T)

#--------- grid arrange shared legend function
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}
