#!/usr/bin/env Rscript

#####after BPnet.R in order to test only the edges of interest
######################GET DATA

args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
BP=args[2]
components=as.numeric(args[3])

library(data.table)
#matrixes in data/ have features in columns
files=list.files("data",full.names=T)
#just want the transcripts of the BP
files=files[grep("transcripts",files,invert=T)]
#what if u had a matrix for each BP in data/???????????????
#maybe if u create the files list instead of reading it
#like paste(subtype,c("CpG","miR","TF",BP))

data=lapply(files,fread)
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],colnames=x$V1))
names(data)=gsub("data/Her2.","",files)

#force order: CpG→miR→TF→BP
data=data[c(1,3,4,2)]####check this with other BPs!!!!!!!!!!!!
names(data)[4]="BP"

######################RANDOMIZED SGCCA
library(mixOmics)
library(igraph)

penalty=read.table(paste(subtype,BP,"descriptors",sep='.'))
oriEdges=fread(paste(subtype,BP,"subnet",sep='.'))
i=paste(oriEdges$source,oriEdges$target)

#Randomizing data
random=lapply(data[1:3],function(x) 
	apply(x,2,function(y) sample(y,length(y),replace=T)))
random$BP=data$BP
#Running SGCCA
temp=wrapper.sgcca(X=random,penalty=c(penalty$CpGs[2],
	penalty$miRNAs[2],penalty$TFs[2],1),scale=T,
	scheme="centroid",ncomp=components)#ncomp to explain 50% of transcripts matrix according to mfa.R
#sum_AVE
#print(rowSums(do.call(rbind,temp$AVE$AVE_X)))
source("function_networkAlt.R")#with no plotting
#Getting network
g=lapply(1:components,function(x) network(temp,comp=list(CpGs=x,
	TFs=x,miRNAs=x,BP=x),blocks=1:4)$gR)#to avoid pdf issues
edges=do.call(rbind,lapply(g,function(x) 
	as.data.frame(cbind(get.edgelist(x),E(x)$weight))))
colnames(edges)=c("source","target","corr")
j=paste(edges$source,edges$target)
print(which(i%in%j))
#aint the same edge coz corr value is not the same

#it's NOT a pval of the corelation, but a
#prob of finding some correlation between the nodes by chance
#pvals=sapply(1:nrow(oriEdges),function(y) 
#	sum(sapply(found,function(x) sum(x==y))>0)/length(found))
#write.table(cbind(oriEdges,pvals),paste(subtype,BP,"subnet",
#	sep='.'),sep='\t',quote=F,row.names=F)