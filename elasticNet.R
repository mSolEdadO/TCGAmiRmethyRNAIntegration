#!/usr/bin/env Rscript
library(caret)
library(methods)
load("/home/msoledad/resultados/EigenScaled.RData")

args = commandArgs(trailingOnly=TRUE)
subtipo=which(names(EigenScaled)==args[1])
subtipo=EigenScaled[[subtipo]]
gen=args[2]
rm(EigenScaled);gc()

coefGrid <-  expand.grid(lambda=10^ seq (3 , -2 , length =10),
			alpha=10^ seq (0, -3 , length =10))
k=5
if(ncol(subtipo)<100){k=3}
trainCtrl <- trainControl("repeatedcv",
			 number = k, #k choose this according to n
			 repeats=round(100/k),#500 for alpha=0.5
			 verboseIter = F,#T if fit fails,
			 allowParallel=T,
			 returnResamp="all")

set.seed(123)
model <- train(y = subtipo[,colnames(subtipo)==gen],
	       x = subtipo[,colnames(subtipo)!=gen],
	       method = "glmnet",
	       trControl = trainCtrl,
	       tuneGrid = coefGrid)

#coefs=as.matrix(coef(model$finalModel, model$bestTune$lambda))
#coefs=as.matrix(coefs[which(coefs>0),])
#colnames(coefs)=paste(model$results$Rsquared[model$results$RMSE==min(model$results$RMSE)],model$results$RMSE[model$results$RMSE==min(model$results$RMSE)],sep="_")
#write.table(coefs,
#	    file=paste(gen,args[1],sep='.'),
#	    quote=F,
#	    sep='\t')
#write.table(model$results,
#	    file=paste(gen,args[1],"results",sep='.'),
#	    quote=F,
#	    sep='\t',
#	    row.names=F)
save(model,file=paste(gen,args[1],"RData",sep='.'))


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
