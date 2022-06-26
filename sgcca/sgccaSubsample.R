#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
iteration=args[2]
#ncomp=as.numeric(args[2])

#library(igraph)
library(mixOmics)
library(data.table)
########################DATA
data=fread(paste(subtype,"eigeNormi",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)

#sample data
i=ncol(data)
data=data[,sample(1:i,round(i/2))]
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")

########################THE SGCCA
penalty=c(CpGs=0.02,transcripts=0.02,miRNAs=0.05)#output of choose_penalty.R
ncomp=nrow(data$miRNAs)-1#the last comp has all loadings>0
final=wrapper.sgcca(X=data,penalty=penalty,scale=F,
	scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
#get selected features
selected=do.call(rbind,lapply(1:3,function(y) 
#y covers the 3 data blocks
	do.call(rbind,lapply(1:ncomp,function(x) 
		cbind("comp"=x,
			"feature"=selectVar(final,comp=x,block=y)[[1]][[1]])))))
#selectVar output is a list of length 2,1st element is a list with name & value)
write.table(selected,paste(subtype,iteration,"selected",sep='.'),sep='\t',
	quote=F,row.names=F)

