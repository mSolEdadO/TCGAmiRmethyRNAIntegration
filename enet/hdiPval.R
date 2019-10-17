#!/usr/bin/env Rscript
library(glmnet)
library(methods)
library(data.table)
library(hdi)
enet=function(y,x,s,ppC,ppE,ppM){
fit <- glmnet(y = y,
	       x = x,
	       alpha = 0.5,
	       lambda = s,
	       standardize=T,
	       penalty.factor=c(rep(ppC,384575),rep(ppE,16475),rep(ppM,433))) 
 m   <- predict(fit, type = "nonzero")
 m[[1]]}


args = commandArgs(trailingOnly=TRUE)
subtipo=args[1]
gen=args[2]
print("Cargo datos")
subtipo=fread(subtipo,sep='\t')
nombres=subtipo$V1
subtipo$V1=NULL
subtipo=t(as.matrix(subtipo))
colnames(subtipo)=nombres
#i=round(nrow(subtipo)*0.8)
#subtipo=subtipo[(i+1):nrow(subtipo),]

#splits=nrow(subtipo)*(nrow(subtipo)-1)/2
fit.multi=multi.split(x=subtipo[,colnames(subtipo)!=gen],
		      y=subtipo[,colnames(subtipo)==gen],
		      model.selector=enet,
		      args.model.selector = list(s=args[3],ppC=args[4],ppE=args[5],ppM=args[6]).
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
