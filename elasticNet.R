#!/usr/bin/env Rscript
library(glmnet)
library(methods)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
subtipo=args[1]
gen=args[2]
print("Cargo datos")
subtipo=fread(subtipo,sep='\t')
nombres=subtipo$V1
subtipo$V1=NULL
subtipo=t(as.matrix(subtipo))
colnames(subtipo)=nombres

model <- glmnet(y = subtipo[,colnames(subtipo)==gen],
	       x = subtipo[,colnames(subtipo)!=gen],
	       alpha = 0.5,
	       lambda = args[3],
	       standardize=T,
	       penalty.factor=c(rep(args[4],384575),rep(args[4],16475),rep(args[6],433))) 

coefs=as.matrix(coef(model))
write.table(coefs[coefs!=0,],
	    file=paste(gen,args[1],"coefs",sep='.'),
	    quote=F,
	    sep='\t')
write.table(paste(args[1],args[2],args[3],args[4],args[5],args[6],model$nulldev,model$dev.ratio,sep=" "),
#nulldev=2*(model with a free parameter per observatio loglikelihood-loglikelihood of the intercept model)
#dev.ratio=fraction of (null) deviance explained=1-dev/nulldev
	    file=paste(gen,args[1],"results",sep='.'),
	    quote=F,
	    sep='\t',
	    row.names=F)
