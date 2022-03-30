#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
net=commandArgs(trailingOnly=TRUE)
library(igraph)
library(tidyverse)

edges=read_tsv(net,col_names=F)
g=graph.data.frame(edges[,1:2])
V(g)$type=substr(V(g)$name,1,1)

#get readable names
methy=read_tsv("../MapMethy.tsv")
methy=methy%>%filter(IlmnID%in%V(g)$name)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
							  'entrezgene_id'),mart)

names=data.frame(cbind("features"=V(g)$name,"name"=V(g)$name))