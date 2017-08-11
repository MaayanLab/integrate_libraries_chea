library(VennDiagram)
library(extrafont)
library(extrafontdb)

area1 = unique(colnames(chea2))
area2 = unique(colnames(encode2))
area3 = unique(colnames(coexp))
area4 = unique(gsub("-.*","",colnames(creedsFisher)))

n12 = intersect(area1,area2)
n13 = intersect(area1,area3)
n14 = intersect(area1,area4)
n23 = intersect(area2,area3)
n24 = intersect(area2,area4)
n34 = intersect(area3,area4)

n123 = intersect(n12,area3)
n124 = intersect(n12,area4)
n134 = intersect(n13,area4)
n234 = intersect(n23,area4)

n1234 = intersect(n123,area4)

grid.newpage()

draw.single.venn(22, category = "Doge People", lty = "solid", 
                 col = "royal blue",fill = "lightpink1", 
                 alpha = 0.5,label.col="cornflower blue",cex=2, 
                 fontfamily = windowsFont("Candara"),cat.col="lightpink1",
                 cat.cex=3,cat.fontfamily = windowsFont("Nyala"))

grid.newpage()
draw.quad.venn(area1 = length(area1), area2 = length(area2), 
               area3 = length(area3), area4 = length(area4), 
               n12 = length(n12), n13 = length(n13), n14 = length(n14), 
               n23 = length(n23), n24 = length(n24), n34 = length(n34),
               n123 = length(n123), n124 = length(n124), n134 = length(n134),
               n234 = length(n234), n1234 = length(n1234),
               category = c("ChEA","ENCODE","Co-expression","CREEDS"), lwd = rep(2,4),
               col = "black",fill = c("dodgerblue","forestgreen",
               "yellow1","red"), alpha = rep(0.5,4),cex = rep(1.5,15),fontfamily = "Arial",
               cat.cex = rep(1.5,4), cat.fontface = rep("bold",4), cat.fontfamily = rep("Arial",4))
