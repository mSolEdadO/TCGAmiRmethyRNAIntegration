#################
#plot intersections
library(igraph)
#matrix subtypes vs function
temp=lapply(enriched,function(x) x%>%distinct(name,subtype)%>%
	table%>%t)#t() so columns are
#count exclusive functions
exclusive=lapply(temp,function(y) 
 rev(table(apply(y[which(rowSums(y==0)==4),],1,paste,collapse=""))))#get all intersections
#$BP
#10000 01000 00100 00010 00001 
#   34     4    18    21    32 
#$KEGG
#10000 01000 00100 00010 00001 
#   22     1     4    12    20 
g=lapply(temp,function(x) crossprod(x))
#$BP
#        subtype
#subtype  Basal Her2 LumA LumB Normal
#  Basal     53    3    9    6     10
#  Her2       3   13    5    1      5
#  LumA       9    5   38    7     11
#  LumB       6    1    7   36      6
#  Normal    10    5   11    6     54
#
#$KEGG
#        subtype
#subtype  Basal Her2 LumA LumB Normal
#  Basal     45    4    7   13     17
#  Her2       4    9    2    2      7
#  LumA       7    2   11    4      6
#  LumB      13    2    4   26      9
#  Normal    17    7    6    9     42
values=lapply(1:2,function(y) lapply(1:5,function(x) 
	c(diag(g[[y]])[x]-exclusive[[y]][x],
		exclusive[[y]][x])))
#graph from the matrix
g1=lapply(g,function(x) 
	graph.adjacency(x,weighted=T,mode="undirected",diag=F))
V(g1$BP)$size=diag(g$BP)
V(g1$KEGG)$size=diag(g$KEGG)
pdf("enrichment.pdf")
lapply(1:2,function(x) 
	plot(g1[[x]],vertex.shape="pie",vertex.pie=values[[x]],
	vertex.size=V(g1[[x]])$size,
	edge.width=E(g1[[x]])$weight*3,
	vertex.pie.color=list(c("gray","red")),
	edge.label=E(g1[[x]])$weight,vertex.frame.color="white",
	vertex.label.cex=1.5,vertex.label.color="black",
	edge.label.cex=1.5,edge.label.color="black",main=names(g1)[x]))
dev.off()
