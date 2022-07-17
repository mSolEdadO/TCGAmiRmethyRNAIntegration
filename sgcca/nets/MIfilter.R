#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
library(tidyverse)
library(igraph)

mi=commandArgs(trailingOnly=TRUE)
edges=read_tsv(mi,col_names=F)

g=graph.data.frame(edges[,1:2],directed=F)
E(g)$MI=edges$X3
g=simplify(g,edge.attr.comb="first")#since MI(a,b)=MI(b,a)
features=V(g)$name
#just gonna ignore CpG-CpG pair?????????????????

#################################real CpG-transcript edges
#library(biomaRt)
methy=read_tsv("../MapMethy.tsv")
methy=methy%>%filter(IlmnID%in%features)
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
#	version=105)
#https://dec2021.archive.ensembl.org
#myannot=getBM(attributes = c("ensembl_gene_id","refseq_ncrna", 
#	"refseq_mrna","mirbase_id"), mart=mart)
#write_tsv(myannot,"myannot")
myannot=read_tsv("myannot")
methy=myannot%>%pivot_longer(-c(1,4),names_to="type",
	values_to="refseq")%>%merge(methy,by="refseq",all.y=T)
methy=methy%>%select(IlmnID,ensembl_gene_id,mirbase_id)%>%
	pivot_longer(-1,values_to="target",names_to="type")%>%
	filter(!is.na(target))
#get MI between every pair. 
#how does it compare with obtained CpG-transcript edges?

#################################real TF-transcript edges
tfs=read_tsv("/home/msoledad/param-rm/TFtargets.tsv")
tfs=tfs%>%separate_rows(target,sep=',',convert=T)%>%
	separate_rows(TF,sep=',',convert=T)%>%filter(target%in%features)

#################################real miR-transcript edges
library(multiMiR)
mirIDs=read_tsv("/home/msoledad/param-rm/miR.ids.map.tsv",skip=1)
#add precursors as possible regulators
temp=mirIDs[,c(1,1)]
colnames(temp)[2]="mature"
temp$mature=gsub("mir","miR",temp$mature)
mirIDs=unique(rbind(mirIDs,temp))
mirIDs=mirIDs[mirIDs$precursor%in%features,]
#get regulatory interactions with miRNAs in feature set
miRtargets=get_multimir(mirna=c(mirIDs$mature,features),
	summary=F,table="validated",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="validated",
	columns=columns(miRtargets),keytype="type")
	colnames(miRtargets)[3]="mature"
miRtargets=merge(miRtargets,mirIDs,by="mature")
#join all pairs
reguEdges=mapply(c, methy[,c("IlmnID","ensembl_gene_id")],
				tfs[,c("TF","target")],
				miRtargets[,c("precursor","target_ensembl")])
gr=graph.data.frame(reguEdges)

###########################################get real MI values
subtype=unlist(strsplit(mi,".",fixed=T))[2]
data=data.table::fread(paste(subtype,"eigeNormi",sep='.'))
data=data[data$V1%in%V(gr)$name,]
