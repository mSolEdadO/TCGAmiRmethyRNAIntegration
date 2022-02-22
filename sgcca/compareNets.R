#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
fun=commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(igraph)
library(ape)
library(ggplot2)
library(RCy3)
#########################LOAD NETS
files=list.files()
files=files[grep(fun,files)]
files=files[grep("tsv",files)]
edges=lapply(files,read_tsv,col_names=F)
names(edges)=sapply(strsplit(files,".",fixed=T),function(x) x[1])

nets=lapply(edges,function(x) {g=graph.data.frame(x[,1:2]);#directed graph
							   E(g)$weight=x$X4;
							   return(g)})
#########################WHICH SUBTYPES ARE MORE ALIKE?
d=sapply(nets,function(x) 
	sapply(nets,function(y) 
		ecount(intersection(x,y))/ecount(graph.union(x,y))))
tree=hclust(as.dist(1-d))
#to build a consensus with all the functions
write.tree(as.phylo(tree),paste(fun,"newick",sep='.'))

#if u wanna plot
png(paste(fun,"tile.png",sep='.'))
cbind(rownames(d),d)%>%as.data.frame%>%
	pivot_longer(-1,names_to="subtype",values_to="ecount")%>%
	ggplot(aes(V1,subtype,fill=as.numeric(ecount)))+geom_tile()+
	scale_fill_viridis_c(name="shared edges")+xlab("")+ylab("")+
	theme(text=element_text(size=18))
dev.off()
#valdría la pena separar por tipo de interacción????
#↑
#↑
#↑
#↑
#########################INTERSECTION CHARACTERIZACION
shared=intersection(nets$Basal,
					nets$Her2,
					nets$LumA,
					nets$LumB,
					nets$Normal)
#DO WEIGHTS CHANGE ACROSS SUBTYPES?
ws=edge_attr(shared)
ws=do.call(rbind,lapply(1:5,function(x) 
	cbind(subtype=names(nets)[x],weight=ws[[x]])))
png(paste(fun,"weights.png",sep='.'))
ws%>%as.data.frame%>%ggplot(aes(x=as.numeric(weight),
	fill=subtype,col=subtype))+geom_density(aes(y=..scaled..),
	alpha=0.5)+xlab("weight")+theme(text=element_text(size=18))+
	scale_fill_viridis_d(option = "plasma")+
	scale_color_viridis_d(option="plasma")
dev.off()
#plot networks
V(shared)$type=substr(V(shared)$name,1,1)
gc=createNetworkFromIgraph(shared,"myIgraph")
setNodeColorMapping(table.column="type",
    mapping.type="discrete",
    colors=c('#5577FF','#FFFFFF','#FF7755'))
setEdgeLineWidthMapping(table.column="type",
    mapping.type="discrete",
    colors=c('#5577FF','#FFFFFF','#FF7755'))
exportImage(filename= paste(fun,"png",sep='.'),type="PNG")
saveSession(filename = paste(fun,"cys",sep='.'))
#WHAT MAKES SPECIAL THIS EDGES?????????
#FISHER EACH TYPE OF REGULATION SHARED VS EXCLUSIVE?

