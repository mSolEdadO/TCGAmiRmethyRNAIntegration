#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)
cells=args[1]
perturbation=args[2]
library(data.table)

###################GET SETS
#library(biomaRt)
#data to translate ensembl ids to symbols
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
#   host="http://nov2020.archive.ensembl.org")
#this version has all the genes, the newest dont
#myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
#             mart=mart)
#get community data
#files=list.files()
#files=files[grep("comm",files)]
#merge all in a data.frame
#communities=lapply(files,read.csv,header=T)
#communities=data.frame(do.call(rbind,lapply(1:4,function(x) 
#   cbind(communities[[x]],paste("clust",x,sep='')))))
#translate to hgnc_symbols
#communities=communities[order(communities$genes),]
#i=table(communities$genes)
#communities$symbol=unlist(sapply(1:length(i),function(x) 
#   rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(i)[x]],
#       i[x])))
#communities$name=paste(communities[,3],"community",communities$community,sep='')
#write.table(communities[,c(1,4,5)],"communities.tsv",sep='\t',quote=F,row.names=F)
communities=read.table("communities.tsv",sep='\t',header=T)
data=fread(paste(cells,"mtrx",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#print(dim(data))
#separate by community
sets=lapply(unique(communities$name),function(x) 
    communities$symbol[communities$name==x])
#list elements must be named for fgsea to work
names(sets)=unique(communities$name)

###################ACTUAL ENRICHMENT
library(fgsea)

gseaPerPertu=function(mtrx,pertu){
i=which(colnames(mtrx)==pertu)  
rank=mtrx[,i]
names(rank)=rownames(mtrx)
rank=rank[order(rank,decreasing=T)]
#print(head(rank));
fgseaRes=fgsea(pathways=sets,stats=rank,nperm=10000);
fgseaRes$leadingEdge=NULL
cbind(pertu,fgseaRes)}
results=gseaPerPertu(data,perturbation)
results=results[results$padj<0.05,]
write.table(results,paste(cells,perturbation,"gsea",sep='.'),sep='\t',quote=F,row.names=F)
