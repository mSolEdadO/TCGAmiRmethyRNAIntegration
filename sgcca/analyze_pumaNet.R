#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
net=commandArgs(trailingOnly=TRUE)
library(igraph)
library(tidyverse)
library(biomaRt)
library(RCy3)

edges=read_tsv(net,col_names=F)
g=graph.data.frame(edges[,1:2])
V(g)$type=substr(V(g)$name,1,1)
V(g)$size=degree(g)
E(g)$width=abs(edges$X4)

#get readable names
methy=read_tsv("../MapMethy.tsv")
methy=methy%>%filter(IlmnID%in%V(g)$name)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
							  'entrezgene_id'),mart=mart)
#CpG names
temp=myannot%>%distinct(hgnc_symbol,entrezgene_id)%>%
			   filter(!is.na(entrezgene_id))%>%
			   merge(methy%>%dplyr::select(entrezgene_id,IlmnID))%>%
			   distinct(hgnc_symbol,IlmnID)
nodes=data.frame("IlmnID"=V(g)$name)
nodes=merge(nodes,temp,by="IlmnID",all.x=T)
#gene names
temp=myannot%>%filter(ensembl_gene_id%in%V(g)$name)%>%
			   dplyr::select(-entrezgene_id)
colnames(nodes)[1]="ensembl_gene_id"
nodes=merge(nodes,temp,by="ensembl_gene_id",all.x=T)
#merge hgnc_symbols for the 2 type of features
colnames(nodes)[1]="feature"
nodes=nodes%>%unite("name",c(hgnc_symbol.x,hgnc_symbol.y),na.rm=T)%>%
	  #paste all the names linked to the same feature
	  group_by(feature)%>%summarise(name=paste(name,collapse=','))
nodes$name[nodes$name==""]=nodes$feature[nodes$name==""]
#order all
if(vcount(g)!=sum(V(g)$name==nodes$feature)){
	stop("names didn't match")
}
nodes=nodes[order(match(nodes$feature,V(g)$name)),]
V(g)$label=nodes$name

gc=createNetworkFromIgraph(g,"myIgraph")
setNodeShapeDefault('ELLIPSE')
setNodeLabelMapping(table.column="label")
setNodeFontSizeDefault(18)
setEdgeLineWidthMapping(table.column="width",
	widths=c(1,20))#pass widths or it'll make huge edges
setNodeColorMapping(table.column="type",
    mapping.type="discrete",
    table.column.values=list('c','E','h'),#automatic mapping fails, so pass this 
    colors=RColorBrewer::brewer.pal(n = 3, name = 'Accent'))
setNodeSizeMapping(table.column="size",
    mapping.type="continuous",
#pass values & sizes or it'll fail    
    table.column.values=range(V(g)$size),
    sizes=c(30,100))
