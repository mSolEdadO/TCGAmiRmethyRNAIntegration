#!/usr/bin/env Rscript

################ libraries #############################
library(data.table)
library(doParallel)
library(hdi)
library(methods)

################ arrange input #########################
#registerDoParallel(cores=3)
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
gen=args[2]
subtype=fread(paste("../",subtype,sep=""),sep='\t')
nombres=subtype$V1
subtype$V1=NULL
#i=round(ncol(subtype)*0.8)
subtype=t(as.matrix(subtype))
colnames(subtype)=nombres
print("done input data")
k=5
if(nrow(subtype)<100){k=3}
registerDoParallel(cores=k)

fit.multi=multi.split(
	y=as.matrix(subtype[,colnames(subtype)==gen]),
	x=as.matrix(subtype[,colnames(subtype)!=gen]),
	model.selector = lasso.cv,
	args.model.selector = list(
				nfolds=k,
				lambda=10^ seq (3, -3 , length=50),
				intercept=F,
				alpha=0.5,
				standardize=T),
	parallel=T,verbose=F,
	ncores=k,
	B=100,ci=F)
write.table(fit.multi$pval.corr,
	    paste(args[1],gen,"pvals",sep='.'),
	    sep='\t',
	    quote=F,
	    col.names=F)

