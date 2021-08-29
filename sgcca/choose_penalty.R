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
	description$nfeat=sapply(resus$loadings,function(x) sum(x!=0))
	description$omic=rownames(description)
	description$penalty=resus$penalty
	colnames(description)[1]="AVE"
return(description)}

results=do.call(rbind,lapply(1:5,function(x) #5 subsamplings're enough???
	describe(data,penalty_cpgs,penalty_transcripts,penalty_mir)))
write.table(results,
	paste(subtype,penalty_cpgs,penalty_transcripts,penalty_mir,
	sep='.'),sep='\t',quote=F,rownames=F)

#############################################
#grid=c(seq(0.01,.1,.01),seq(0.3,.9,.2))
#sub=readLines("condor.seed")
results=pblapply(grid,function(x) sapply(grid,function(y) sapply(grid,function(z) describe(data,x,y,z))))