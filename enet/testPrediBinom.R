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
	       intercept=F,
	       penalty.factor=c(rep(ppC,384575),rep(ppE,16475),rep(ppM,433))) 
 predict(fit, type = "coefficients")[,1]}


args = commandArgs(trailingOnly=TRUE)
subtipo=args[1]
gen=args[2]
subtipo=fread(subtipo,sep='\t')
nombres=subtipo$V1
subtipo$V1=NULL
subtipo=t(as.matrix(subtipo))
colnames(subtipo)=nombres
#i=round(nrow(subtipo)*0.8)
#subtipo=subtipo[(i+1):nrow(subtipo),]
k=5
if(nrow(subtipo)<100){k=3}

stability=do.call(cbind,lapply(1:1000,function(x) {
 i=sample(1:nrow(subtipo),round(nrow(subtipo)/k));
 print(x);
 enet(x=subtipo[i,colnames(subtipo)!=gen],
     y=subtipo[i,colnames(subtipo)==gen],
     s=args[3],
     ppC=args[4],
     ppE=args[5],
     ppM=args[6])}))
stability=rowSums(stability)

write.table(stability,
	    paste(args[1],gen,"pvalSol",sep='.'),
	    sep='\t',
	    quote=F,
	    col.names=F)

