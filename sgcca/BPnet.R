#!/usr/bin/env Rscript
library(data.table)
library(igraph)
library(biomaRt)

BP=commandArgs(trailingOnly=TRUE)
edges=fread(paste(BP,"net",sep='.'))

########################focus on BP nodes##################
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
print(paste("Checking",BP,"in biomaRt"))#this takes a loooong time
geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),
	filters = 'go',values = BP, mart = mart)
#using BP as filter outputs other GOs eitherway
geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP])
#verify BP nodes are in col2
print(sum(edges$source%in%geneSet))
print(sum(edges$target%in%geneSet))
edges=edges[edges$target%in%geneSet,]

########################collapse CpGs per CGI##################
methy=fread("jhu-usc.edu_BRCA.HumanMethylation450.6.lvl-3.TCGA-AO-A03N-01B-11D-A10N-05.gdc
_hg38.txt")
colnames(methy)[1]="ID"

#verify CpGs are in col1
print(table(substr(edges$source,1,1)))
print(table(substr(edges$target,1,1)))
#separate edges with & without CpGs
methy=methy[methy$ID%in%edges$source,]
temp=edges[substr(edges$source,1,1)=="c",]
edges=edges[substr(edges$source,1,1)!="c",]
#add CGI info
temp=temp[order(temp$source),]
i=table(temp$source)
temp$CGI=unlist(sapply(1:length(i),function(x) 
	rep(methy$CGI_Coordinate[methy$ID==names(i)[x]],i[x])))

g=graph.data.frame(temp[,c(2,4)],directed=F)#targets to CGI
E(g)$cor=as.numeric(temp$corr)
#magic happens here
g1=simplify(g,edge.attr.comb=list(cor=function(x) 
	median(abs(x)),name="ignore"))
print(paste("Components in subnetwork:",components(g1)$no),sep=' ')
edgesAlt=as.data.frame(get.edgelist(g1))
edgesAlt$corr=E(g1)$cor
edgesAlt=rbind(as.matrix(edgesAlt),as.matrix(edges))