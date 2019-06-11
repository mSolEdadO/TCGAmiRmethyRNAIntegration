#!/usr/bin/env Rscript
library(caret)
library(methods)
load("../resultados/Eigenscaled.RData")#= MFA output: concatenated matrix normalized to weight the omics per variance not size

args = commandArgs(trailingOnly=TRUE)
subtipo=which(names(Eigenscaled)==args[1])
subtipo=Eigenscaled[[subtipo]] #fit this subtype matrix
gen=args[2]#to this PAM50 gene expression

#over this parameters
coefGrid <-  expand.grid(lambda=10^ seq (2 , -2 , length =10),
			alpha=10^ seq (0, -3 , length =10))
#with this training
trainCtrl <- trainControl("repeatedcv",
			 number = 10, #k choose this according to n
			 repeats=100,#?????
			 verboseIter = F,#T if fit fails,
			 allowParallel=T,
			 returnResamp="all")

set.seed(123)
model <- train(y = subtipo[,colnames(subtipo)==gen],
	       x = subtipo[,colnames(subtipo)!=gen],
	       method = "glmnet",
	       trControl = trainCtrl,
	       tuneGrid = coefGrid)

#save coefficients and model statistics
coefs=as.matrix(coef(model$finalModel, model$bestTune$lambda))#sif
colnames(coefs)=paste(model$bestTune,collapse='_')
write.table(coefs,
	    file=paste(gen,args[1],sep='.'),
	    quote=F,
	    sep='\t')
write.table(model$results,
	    file=paste(gen,args[1],"results",sep='.'),
	    quote=F,
	    sep='\t')
write.table(model$resample,
	    file=paste(gen,args[1],"resample",sep='.'),
	    quote=F,
	    sep='\t')

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
