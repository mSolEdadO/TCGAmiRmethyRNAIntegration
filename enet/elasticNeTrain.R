#!/usr/bin/env Rscript
library(caret)
library(methods)
library(data.table)
library(doParallel)
#load("/home/msoledad/resultados/EigenScaled.RData")
registerDoParallel(cores=3)

args = commandArgs(trailingOnly=TRUE)
#subtipo=which(names(EigenScaled)==args[1])
#subtipo=EigenScaled[[subtipo]]
subtipo=args[1]
gen=args[2]
print("Cargo datos")
subtipo=fread(paste("../",subtipo,sep=""),sep='\t')
nombres=subtipo$V1
subtipo$V1=NULL
i=round(ncol(subtipo)*0.8)
subtipo=t(as.matrix(subtipo)[,1:i])#training subset should not be used for testing
colnames(subtipo)=nombres
#rm(EigenScaled);gc()

coefGrid <-  expand.grid(lambda=10^ seq (3 , -2 , length =100),
			alpha=0.5)
k=5
if(nrow(subtipo)<100){k=3}
trainCtrl <- trainControl("repeatedcv",
			 number = k, #k choose this according to n
			 repeats=600/k,#200 for the small training set, 
			 verboseIter = T,#T if fit fails,
			 allowParallel=T,
			 returnResamp="all")

set.seed(123)
model <- train(y = subtipo[,colnames(subtipo)==gen],
	       x = subtipo[,colnames(subtipo)!=gen],
	       method = "glmnet",
	       trControl = trainCtrl,
	       tuneGrid = coefGrid,
	       standardize=T,
	       penalty.factor=c(rep(0.1,384575),rep(1,16475),rep(0.5,433)))#variate

#coefs=as.matrix(coef(model$finalModel, model$bestTune$lambda))
#coefs=as.matrix(coefs[which(coefs>0),])
#colnames(coefs)=paste(model$results$Rsquared[model$results$RMSE==min(model$results$RMSE)],model$results$RMSE[model$results$RMSE==min(model$results$RMSE)],sep="_")
#write.table(coefs,
#	    file=paste(gen,args[1],"coefs",sep='.'),
#	    quote=F,
#	    sep='\t')
write.table(model$results,
	    file=paste(gen,args[1],"results",sep='.'),
	    quote=F,
	    sep='\t',
	    row.names=F)
#save(model,file=paste(gen,args[1],"RData",sep='.'))


#plus elasticNet & elasticNet.sub -> paralell
########################################
#elasticNet should have this
#!/bin/bash

## run R, with the name of your  R script
#../R/bin/Rscript elasticNet.R $1 $2
########################################
#elasticNet.sub should have this
# Opciones generales de HTCondor.
#universe = vanilla
#
#initialdir = /castle/msoledad/parallel-caret
#should_transfer_files = NO
#getenv = True
#
# Recursos necesarios para ejecutar el trabajo.
#request_cpus = 1
#request_memory = 3GB
#request_disk = 1GB

#per gene & subtype
#executable = elasticNet
#arguments =  LumA ENSG00000141738
#log =LumA.ENSG00000141738.condor_log
#error =LumA.ENSG00000141738.condor_err
#queue
#################################
#sub=readLines("temp")
#listos=list.files()
#listos=listos[grep("results",listos)]
#listos=sapply(strsplit(listos,".",fixed=T),function(x) paste(x[2],x[1]))
#sub=sub[!1:length(sub)%in%sapply(listos,function(x) grep(x,sub))]
#writeLines(sub,"temp")
