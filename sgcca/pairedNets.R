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
if(length(grep("hsa",fun))>0){
 funData=read_tsv("../KEGG-allFeatures.enrichment")
 library(biomaRt)
 mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
 myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
							  'entrezgene_id'),mart=mart)
 funData=funData%>%separate_rows(geneID,sep='/')
 colnames(myannot)[3]="geneID"
 funData=merge(funData,myannot,by="geneID")
 features=funData%>%filter(ID%in%fun)%>%group_by(subtype)%>%
 			group_map(~unique(.x$hgnc_symbol))
}else{
 funData=read_tsv("../BP-allFeatures.enrichment")
 features=funData%>%filter(ID%in%fun)%>%group_by(subtype)%>%
	group_map(~unique(unlist(strsplit(unlist(.x$geneID),"/"))))
}
names(features)=funData%>%filter(ID%in%fun)%>%group_by(subtype)%>%
				group_keys%>%unlist
#ego_graph for eniched features
featID=lapply(names(nets),function(y) 
	unique(unlist(lapply(features[[y]],function(x) 
		#get vid with grep coz CpGs labels are comma separated
		grep(x,V(nets[[y]])$label)))))
names(featID)=names(nets)
subn=lapply(names(nets),function(x) 
	subgraph.edges(nets[[x]],
		E(nets[[x]])[from(featID[[x]])|to(featID[[x]])]))

#paired comparissons
gc=createNetworkFromIgraph(subn[[1]],
	title=paste(names(nets)[1],paste(fun,collapse=',')))
gd=createNetworkFromIgraph(subn[[2]],
	title=paste(names(nets)[2],paste(fun,collapse=',')))
#merge them manually in cytoscape coz igraph duplicates attributes
#importVisualStyles("mystyle.xml")
#[1] "default_0" apply this style 
#take care of intersection nodes
saveSession()

####################check edges on dbs
#get gene description from biomaRt for all nodes
library(biomaRt)
names=data.frame(unique(do.call(rbind,sapply(subn,function(x)
 cbind("name"=V(x)$name,"label"=V(x)$label)))))
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
#https://dec2021.archive.ensembl.org
myannot=getBM(attributes = c("hgnc_symbol","wikigene_description"),
	filter="hgnc_symbol",values=names$label,mart=mart)

#rentrez per node and function
library(rentrez)
#search only the nodes that aren't responsible for the enrichment
query=names$label[!names$label%in%unlist(features)]
query1=funData%>%filter(funData$ID%in%fun)%>%
 dplyr::select(Description)%>%unique%>%list
knownfun=lapply(query,function(x) lapply(query1,function(y) 
	entrez_search(db="pubmed",term=paste(x,y,sep=' AND '))))
lapply(knownfun,function(x) lapply(x,function(y) y$ids[y$count>0]))

edges=data.frame(unique(do.call(rbind,lapply(subn,function(x) 
	get.edgelist(x)))))
colnames(edges)[1]="name"
edges=merge(edges,names,by="name")[,2:3]
colnames(edges)[2]="pair1"
colnames(edges)[1]="name"
edges=unique(merge(edges,names,by="name")[,2:3])
colnames(edges)[2]="pair2"
known=apply(edges,1,function(x) 
	entrez_search(db="pubmed",
				  term=paste(x[1],x[2],sep=' AND ')))
known=known[sapply(known,function(x) x$count)>0]
sapply(known,function(x) x$ids)

#multimiR miRNA-target prediction
library(multiMiR)
query=query[grep("mir",query,ignore.case=T)]
query=paste("hsa",gsub("miR","mir",query),sep='-')
mirIDs=read_tsv("/home/msoledad/param-rm/miR.ids.map.tsv",skip=1)
mirIDs=mirIDs[mirIDs$precursor%in%query,]
miRtargets=get_multimir(mirna=c(query,mirIDs$mature),
	summary=F,table="predicted",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="predicted",
	columns=columns(miRtargets),keytype="type")
miRtargets[miRtargets$target_symbol%in%names$label,]
#what about TFs????

#########################WHICH SUBTYPES ARE MORE ALIKE?
#########################WHICH TYPE OF EDGES ARE MORE SHARED
#########################EXCLUSIVE EDGES

#WHAT MAKES SPECIAL THIS EDGES?????????



