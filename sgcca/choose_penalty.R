#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
penalty_cpgs=args[2]
penalty_transcripts=args[3]
penalty_mir=args[4]

######################DATA TO LIST OF MATRIXES PER MOLECULAR LEVEL
library(data.table)
data=fread(paste(subtype,"normalized",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
n=ncol(data)
subn=round(n*.5)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")

######################CHECK SPARSITY VALUES
#sparsity parameters are chosen for each of the 10 MCCV iterations
#using an internal 5-fold CV loop: the parameters that minimize the
#prediction error [Tenenhaus2014]

library(mixOmics)

#take model descriptors<-----------------recycled
describe=function(data,pc,pt,pm){
	#subsample observations
	data=lapply(data,function(x) x[sample(1:n,subn),])
	resus=wrapper.sgcca(data,penalty=c(pc,pt,pm),scale=T,
		scheme="centroid")
	#get results description
	description=as.data.frame(do.call(rbind,resus$AVE$AVE_X))
	description$nfeatures=sapply(resus$loadings,function(x) sum(x!=0))
	description$omic=rownames(description)
	description$penalty=resus$penalty
	colnames(description)[1]="AVE"
return(description)}

#this block is for cluster
#results=do.call(rbind,lapply(1:10,function(x) #are 10 subsamplings enough???
#	describe(data,penalty_cpgs,penalty_transcripts,penalty_mir)))
#write.table(results,
#	paste(subtype,penalty_cpgs,penalty_transcripts,penalty_mir,
#	sep='.'),sep='\t',quote=F,rownames=F)
#use this for making the submit file
#grid=seq(0.01,0.9,length=10)#first values tested
#sub=readLines("condor.seed")
#args=as.character(sapply(grid,function(x) 
#	sapply(grid,function(y) sapply(grid,function(z)
#	paste("Her2",x,y,z,sep=' ')))))
#temp=lapply(args,function(x) gsub("Her2 0.01 0.01 0.01",x,sub))
#writeLines(unlist(temp),"temp")
#cat temp|perl -pe 's/\t/\n/g'>temp1
#cat condor.header temp1>penalty.sub

#this block is for your lap
grid=seq(0.01,0.1,0.01)
resus=pblapply(grid,function(x) lapply(grid,function(y) 
	lapply(grid,function(z) describe(data,x,y,z))))
resus=do.call(rbind,lapply(resus,function(x) 
	do.call(rbind,lapply(x,function(y) do.call(rbind,y)))))
write.table("penalty_search.tsv",sep='\t',quote=F,row.names=F)


