#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
BP=args[2]
components=as.numeric(args[3])
print(components)

######################FORMAT DATA
print("Building data")
library(data.table)
data=fread(paste(subtype,"normalized",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
tfs=read.table("TFtargets.tsv",header=T,sep='\t')#https://github.com/slowkow/tftargets
tfs=unique(unlist(strsplit(as.character(tfs$TF),',')))
data$TFs=data$transcripts[,colnames(data$transcripts)%in%tfs]

#keep only the transcripts for a biological process
library(biomaRt)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
print(paste("Checking",BP,"in biomaRt"))#this takes a loooong time
geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),filters = 'go',
 values = BP, mart = mart)
#using BP as filter outputs other GOs eitherway
geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP])
data$BP=data$transcripts[,colnames(data$transcripts)%in%geneSet]
data$transcripts=NULL
#avoid selecting TFs that are transcripts from BP
data$TFs=data$TFs[,!colnames(data$TFs)%in%colnames(data$BP)]

######################RANDOMIZED SGCCA
library(mixOmics)
library(igraph)

random_sgcca=function(){
	print("Randomizing data")
	random=lapply(data[1:3],function(x) 
		apply(x,2,function(y) sample(y,length(y),replace=T)))
	random$BP=data$BP
	print("Running SGCCA")
	temp=wrapper.sgcca(X=random,penalty=c(penalty$CpGs[2],
		penalty$miRNAs[2],penalty$TFs[2],1),scale=T,
		scheme="centroid",ncomp=components)#ncomp to explain 50% of transcripts matrix according to mfa.R
	print("sum_AVE")
	print(rowSums(do.call(rbind,temp$AVE$AVE_X)))
	source("function_networkAlt.R")#with no plotting
	print("Getting network")
	g=lapply(1:components,function(x) network(temp,comp=list(CpGs=x,
		TFs=x,miRNAs=x,BP=x),blocks=1:4)$gR)#to avoid pdf issues
	edges=do.call(rbind,lapply(g,function(x) 
		as.data.frame(cbind(get.edgelist(x),E(x)$weight))))
	colnames(edges)=c("source","target","corr")
	j=paste(edges$source,edges$target)
return(which(i%in%j))}#aint the same edge coz corr value is not the same

penalty=read.table(paste(subtype,BP,"descriptors",sep='.'))
oriEdges=fread(paste(subtype,BP,"BP_net",sep='.'))
i=paste(oriEdges$source,oriEdges$target)
#how many times should I apply it??????
found=lapply(1:10,function(x) random_sgcca())
#it's NOT a pval of the corelation, but a
#prob of finding some correlation between the nodes by chance
pvals=sapply(1:nrow(oriEdges),function(y) 
	sum(sapply(found,function(x) sum(x==y))>0)/length(found))
write.table(cbind(oriEdges,pvals),paste(subtype,BP,"descriptors",
	sep='.'),sep='\t',quote=F,row.names=F)