#!/usr/bin/env Rscript

library(mixOmics)
library(data.table)

#sparsity parameters are chosen for each of the 10 MCCV iterations
#using an internal 5-fold CV loop: the parameters that minimize the
#prediction error [Tenenhaus2014]
params_searcher=function(subtype,lCpG,ltrnscri,lmir,count){
	size=nrow(subtype$miRNAs)
	i=sample(1:size,round(size/2))
	#test in half the available data
	data=lapply(subtype,function(y) y[i,])
	#ncomp=1 for training
	results=wrapper.sgcca(data,penalty=c(lCpG,ltrnscri,lmir),scale=T)
	#save only the model's size
	evar=as.data.frame(do.call(rbind,results$explained_variance))
	features=as.data.frame(do.call(rbind,lapply(results$loadings,
		function(x) sum(x!=0))))#number of selected features
	out=cbind(evar,features,results$penalty)
	out$omic=rownames(out)
	colnames(out)[1:3]=c("explained_variance","nfeatures","sparsity")
return(out)}

args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
lCpG=as.numeric(args[2])
ltrnscri=as.numeric(args[3])
lmir=as.numeric(args[4])
#load data
data=fread(paste(subtype,"normalized",sep='.'))
#op1:cca=wrapper.sgcca(raw,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op2:cca.n=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=F)
#op3:cca.n1=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op1 & op3 give the same results BUT op1 ends (cca$crit) faster
#is raw data better???what if such covergence depends on penalties???
#selectVar() is the same than loadings != 0, NOT loadings.star != 0
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")

#grid=seq(.01,.9,length=20)
#test several combinations
#temp=do.call(rbind,lapply(grid,function(x) 
#	do.call(rbind,lapply(grid,function(y) 
#        do.call(rbind,lapply(grid,function(z) 
#        	{do.call(rbind,lapply(1:10,function(i) 
#                params_searcher(data,x,y,z)));
#        	print(paste(x,y,z,sep=','))}))))))#how far am I?
#test 1 set of sparsity values
combi=do.call(rbind,lapply(1:20,function(i) 
	params_searcher(data,lCpG,ltrnscri,lmir)))
write.table(combi,paste(subtype,lCpG,ltrnscri,lmir,"params",sep='_'),
	sep='\t',quote=F,row.names=F)
####run it for as many combinations u wanna

#u can build condor submit files like:
#sub=readLines("temp")#.sub with field separated by \t
#sub=unlist(strsplit(sub,"\t"))
#grid=seq(.01,.9,length=10)
#combination of sparsity values
#args=as.character(sapply(grid,function(x) 
#	sapply(grid,function(y) 
#		sapply(grid,function(z) paste(x,y,z,sep=' ')))))
#.sub line per combination
#args=sapply(1:length(args),function(i)
#	gsub("0.1 0.1 0.1",args[i],sub[2]))
#sub=sapply(1:length(args),function(x) 
#	paste(sub[1],args[x],sub[3],sub[4],sub[5],sub[6],sep='_'))
#writeLines(sub,"temp")#add the lines on the needed resources

