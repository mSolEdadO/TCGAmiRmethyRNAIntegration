#!/usr/bin/env Rscript
library(glmnet)
library(methods)
library(data.table)
library(hdi)
library(doParallel)
#load("/home/msoledad/resultados/EigenScaled.RData")
registerDoParallel(cores=3)

args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
gen=args[2]
print("Cargo datos")
subtype=fread(subtype,sep='\t')
nombres=subtype$V1
subtype$V1=NULL
subtype=t(as.matrix(subtype))
colnames(subtype)=nombres
i=round(nrow(subtype)*0.8)
training=subtype[1:i,]

k=5
if(nrow(training)<100){k=3}
fit.multi=multi.split(x=training[,colnames(training)!=gen],
		      y=training[,colnames(training)==gen],
		      #model.selector=enet,
		      args.model.selector = list(nfolds=k,
		      							 intercept=F,
		      							 alpha=0.5,
		      							 standardize = TRUE,
		      							 penalty.factor=c(rep(args[4],384575),rep(args[5],16475),rep(args[6],433)),
		      							 lambda=10^ seq (3 , -2 , length =50)),#length is usually 100
		      verbose=T,
		      #return.selmodels=T,
		      #return.nonaggr=T,
		      ci=F,
		      parallel=T,
		      ncores=4,
		      B=50
		      )

write.table(fit.multi$pval.corr,
	    paste(args[1],gen,"pvals",sep='.'),
	    sep='\t',
	    quote=F,
	    col.names=F)
