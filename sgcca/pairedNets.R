#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
net=args[1]
fun=args[2]

library(tidyverse)
library(igraph)
library(RCy3)

#########################LOAD NETS
files=list.files()
files=files[grep(net,files)]
files=files[grep("cys",files)]
nets=lapply(files,function(x) {openSession(x);
							   createIgraphFromNetwork()})
#exportVisualStyles("mystyle.xml")
names(nets)=sapply(strsplit(files,".",fixed=T),function(x) x[2])
################GET NODES OF INTEREST
if(length(grep("GO",fun))>0){
 funData=read_tsv("../BP-allFeatures.enrichment")
}else{
 funData=read_tsv("../KEGG-allFeatures.enrichment")
}
features=funData%>%filter(ID==fun)%>%group_by(subtype)%>%
	group_map(~unique(unlist(strsplit(unlist(.x$geneID),"/"))))
names(features)=funData%>%filter(ID==net)%>%group_by(subtype)%>%
				group_keys%>%unlist
#ego_graph for eniched features
features=lapply(names(nets),function(y) 
	unique(unlist(lapply(features[[y]],function(x) 
		#get vid with grep coz CpGs labels are comma separated
		grep(x,V(nets[[y]])$label)))))
names(features)=names(nets)
subn=lapply(names(nets),function(x) 
	subgraph.edges(nets[[x]],
		E(nets[[x]])[from(features[[x]])|to(features[[x]])]))

#paired comparissons
gc=createNetworkFromIgraph(subn[[1]],title=paste(names(nets)[1],fun))
gd=createNetworkFromIgraph(subn[[2]],title=paste(names(nets)[2],fun))
#merge them manually in cytoscape coz igraph duplicates attributes
#importVisualStyles("mystyle.xml")
#[1] "default_0" apply this style 
#take care of intersection nodes
saveSession()

#########################WHICH SUBTYPES ARE MORE ALIKE?
#########################WHICH TYPE OF EDGES ARE MORE SHARED
#########################EXCLUSIVE EDGES

#WHAT MAKES SPECIAL THIS EDGES?????????



