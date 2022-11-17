#!/usr/bin/env Rscript

library(tidyverse)
library(biomaRt)

#get data
files=list.files()
files=files[grep("stable",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".stable","",files)

#map every omic id to entrez
methy=read_tsv("../MapMethy.tsv")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=105)
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"mirbase_id"), mart=mart)
#prepare to merge with sets
myannot=myannot%>%pivot_longer(-entrezgene_id,names_to="DB",
	values_to="variable")
myannot=myannot%>%dplyr::select(entrezgene_id,variable)%>%
	filter(variable!="")
methy=methy%>%dplyr::select(entrezgene_id,IlmnID)
colnames(methy)[2]="variable"
myannot=rbind(myannot,methy)
#add entrezgene_id to selected features
sets=lapply(sets,function(x) merge(x,myannot,all.x=T,by="variable"))
#all the transcripts serve for enrichment, some miRNAs & CpGs don't
#sapply(setsAlt,function(x) table(x$omic[is.na(x$entrezgene_id)]))
#       Basal Her2  LumA  LumB Normal
#CpGs   14123 4038 41035 13250  10200
#miRNAs     6    1    25     7      3

#########################ACTUAL ENRICHMENT#########################
library(clusterProfiler)
library(org.Hs.eg.db)

complexset=do.call(rbind,lapply(1:5,function(x) 
	cbind(names(sets)[x],sets[[x]][,c("component","entrezgene_id","omic")])))
colnames(complexset)[1]="subtype"
complexset=complexset%>%filter(entrezgene_id!="NA")
complexset$entrezgene_id=as.character(complexset$entrezgene_id)

BPenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichGO",
	OrgDb=org.Hs.eg.db,
	ont="BP",
	readable=T,
	pAdjustMethod = "fdr",
    pvalueCutoff  = 0.01)#slooow
write_tsv(as.data.frame(BPenrich),"BP.enrichment")

KEGGenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
write_tsv(as.data.frame(KEGGenrich),"KEGG.enrichment")                  
