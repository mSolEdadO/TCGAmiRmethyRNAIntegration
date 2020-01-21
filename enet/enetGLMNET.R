#!/usr/bin/env Rscript

################ libraries #############################
library(data.table)
library(glmnet)
library(doParallel)
library(methods)
library(caret)

################ arrange input #########################
#registerDoParallel(cores=3)
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
gen=args[2]
subtype=fread(paste("../",subtype,sep=""),sep='\t')
nombres=subtype$V1
subtype$V1=NULL
i=round(ncol(subtype)*0.8)
subtype=t(as.matrix(subtype))
colnames(subtype)=nombres
training=subtype[1:i,]
testing=subtype[(i+1):nrow(subtype),]
k=5
if(nrow(training)<100){k=3}
registerDoParallel(cores=k)
print("input done")

################ model fit #############################
cvfit=function(gen,k){
print(dim(training[,colnames(training)!=gen]))
print(length(training[,colnames(training)==gen]))
cv.glmnet(
	y =training[,colnames(training)==gen],
	x =training[,colnames(training)!=gen],
	lambda=10^ seq (3, -3 , length=50),
	nfolds=k,
	parallel=T,
	intercept=F,
	alpha=0.5,
	standardize=T,
	weights=rep(1,nrow(training)))}

#base fit
print("start main fit")
#seed(123)
model=cvfit(gen,k)

#choose the lambda with the smallest MSE
print("start repeated cv")
cvm=sapply(1:99,function(i) cvfit(gen,k)$cvm)
cvm=cbind(cvm,model$cvm)
bestLambda=model$lambda[which.min(rowMeans(cvm))]

################ best lambda predictors ################
coefs=coef(model,s=bestLambda)
coefs1=matrix(coefs[which(coefs!=0)],ncol=1)
rownames(coefs1)=rownames(coefs)[which(coefs!=0)]
ajuste=RMSE(predict(model,testing[,colnames(testing)!=gen]),testing[,colnames(testing)==gen])
colnames(coefs1)=paste("lambda",bestLambda,"RMSE",ajuste,sep=':')
write.table(coefs1,
	    file=paste(gen,args[1],"coefs",sep='.'),
	    quote=F,
	    sep='\t')
