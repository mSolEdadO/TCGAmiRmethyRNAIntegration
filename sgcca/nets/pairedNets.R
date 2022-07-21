#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
net=commandArgs(trailingOnly=TRUE)
id=unlist(strsplit(net,'.',fixed=T))[1]
subty=unlist(strsplit(net,'.',fixed=T))[2]

library(tidyverse)
library(igraph)
library(RCy3)

#########################LOAD NETS
openSession(net)
g=createIgraphFromNetwork()

#################################################HIGHLIGHT TFS
tfs=read_tsv("/home/msoledad/param-rm/TFtargets.tsv",show_col_types=F)
tfs=tfs%>%separate_rows(TF,sep=',',convert=T)%>%
	filter(TF%in%V(g)$name)%>%distinct(TF)%>%unlist
sapply(tfs,function(x) {
	setNodeBorderColorBypass(x, '#FF0000');
	setNodeBorderWidthBypass(x, 8)})

#######################################3HIGHLIGHT FUNCTIONAL NODES
#get the cluster that contains the function
chosen=read_tsv("../subsamples/chosen.tsv")
chosen=chosen%>%filter(ID==id&subtype==subty)
cl=chosen$group
funs=chosen$Description

#get related functions
groups=read_tsv("../subsamples/Groups_per_component.tsv")
groups=groups%>%filter(subtype==subty,group==cl)
if(nrow(groups)>1){
	funs=groups%>%select(Description)%>%unlist}

#get "functional" nodes
enriched=list(BP=read_tsv("../subsamples/BP.enrichment"),
	KEGG=read_tsv("../subsamples/KEGG.enrichment"))
known_genes=lapply(enriched,function(x) 
	x%>%filter(Description%in%funs&subtype==subty)%>%
	dplyr::select(Description,geneID)%>%separate_rows(geneID,sep='/'))
#get symbols for KEGG
if(nrow(known_genes$KEGG)>0){
 library(biomaRt)
 mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
 myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
							  'entrezgene_id'),mart=mart)
 colnames(myannot)[3]="geneID"
 known_genes$KEGG=merge(known_genes$KEGG,myannot,by="geneID",all.x=T)%>%
 	dplyr::select(Description,hgnc_symbol)
 colnames(known_genes$KEGG)[2]="geneID"
 known_genes=do.call(rbind,known_genes)
}else{
 known_genes=known_genes$BP	
}

known_genes=known_genes%>%filter(geneID%in%unique(unlist(strsplit(V(g)$label,","))))
gk=createNetworkFromIgraph(graph.data.frame(known_genes),"known")
mergeNetworks(c(subty,"known"),"merged",nodeKeys=c("label","name"))
#fix names
saveSession()



####################check edges on dbs
#get gene description from biomaRt for all nodes
library(biomaRt)
names=data.frame(cbind("name"=V(g)$name,"label"=V(g)$label))
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
#https://dec2021.archive.ensembl.org
myannot=getBM(attributes = c("hgnc_symbol","wikigene_description"),
	filter="hgnc_symbol",values=names$label,mart=mart)

#rentrez per node and function
library(rentrez)

#m=subgraph.edges(g,E(g)[from(names$name[names$label%in%known_genes$geneID])|to(names$name[names$label%in%known_genes$geneID])])
#query=names%>%filter(name%in%V(m)$name&!label%in%known_genes$geneID)%>%dplyr::select(label)%>%unlist
#search only the nodes that aren't responsible for the enrichment
query=names$label[!names$label%in%known_genes$geneID]
query=unique(unlist(strsplit(query,",")))
knownfun=lapply(funs,function(x) lapply(query,function(y) 
	entrez_search(db="pubmed",term=paste(x,y,sep=' AND '))))
lapply(knownfun,function(x) lapply(x,function(y) y$ids[y$count>0]))

edges=data.frame(get.edgelist(g))
colnames(edges)[2]="name"
edges=merge(edges,names,by="name")[,2:3]
colnames(edges)[2]="pair1"
colnames(edges)[1]="name"
edges=unique(merge(edges,names,by="name")[,2:3])
colnames(edges)[2]="pair2"
edges=edges%>%separate_rows(pair2,sep=',',convert=T)%>%separate_rows(pair2,sep=',',convert=T)
known=apply(edges,1,function(x) 
	entrez_search(db="pubmed",
				  term=paste(x[1],x[2],sep=' AND ')))
known=known[sapply(known,function(x) x$count)>0]
sapply(known,function(x) x$ids)

#multimiR miRNA-target prediction
library(multiMiR)
query=query[grep("mir",query,ignore.case=T)]
if(length(query)>0){
query=paste("hsa",gsub("miR","mir",query),sep='-')
mirIDs=read_tsv("/home/msoledad/param-rm/miR.ids.map.tsv",skip=1)
mirIDs=mirIDs[mirIDs$precursor%in%query,]
miRtargets=get_multimir(mirna=c(query,mirIDs$mature),
	summary=F,table="predicted",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="predicted",
	columns=columns(miRtargets),keytype="type")
miRtargets[miRtargets$target_symbol%in%names$label,]
}
#what about TFs????

########################WHICH TYPE OF EDGES ARE MORE SHARED
#########################EXCLUSIVE EDGES

#down 1D91C0
#up E31A1C

