###########################libraries###############################################
library(data.table)
library(caret)
library(doParallel)
#registerDoParallel(cores=2)
#ridge regression function
ridge=function(expression,inhibition){
	#parameters to fit: lambda=penalization chosen by cv, alpha=0 since this is a ridge regression
	coefGrid <-  expand.grid(lambda=seq (0,50, length =100),#length is usually 100
			alpha=0)
	set.seed(123)
	#cv parameters
	trainCtrl <- trainControl("repeatedcv",
			 number = 3, #min number of folds, 3 since n is small 
			 repeats=100,#the more the better li
			 #verboseIter = T,#T if fit fails,
			 allowParallel=T)
	set.seed(567)
	#actual fit
	model <- train(y = inhibition,
	       x = expression,
	       method = "glmnet",
	       trControl = trainCtrl,
	       tuneGrid = coefGrid,
	       standardize=T,
	       intercept=F)
	return(model)}

###########################data###############################################
#perturbations comparable across cell lines
X=fread("/labs/csbig/DREAMCHALLENGE_CTD2/results/z_matrix.tsv")
#only these genes matter
targets=read.csv("/labs/csbig/DREAMCHALLENGE_CTD2/raw/TARGETS/target_list.csv")
#predicted Y
effects=fread("/labs/csbig/DREAMCHALLENGE_CTD2/results/predicted_effect.tsv")
#concentrations=fread("/labs/csbig/DREAMCHALLENGE_CTD2/results/experimentAnnotation_df.tsv")

###########################formating data#####################################
names=X$name
X=as.matrix(X[,2:ncol(X)])
rownames(X)=names
#only treatments matter?
X=X[,grep("(untreated|dmso)",colnames(X),perl=T,invert=T)]
#perturbation's order should match response order
X=t(X[,order(sapply(strsplit(colnames(X),"_"),function(x) x[2]))])
#treatments are
#sum(substr(rownames(X)[seq(1,704,22)],1,7)%in%effects$tx) 
#cell lines also are in order
#sum(sapply(strsplit(rownames(X)[seq(1,nrow(X),2)],"_"),function(x) x[5])%in%effects$cell_line)==nrow(effects)
#get a matrix per tx
#only kept drugable genes
X=X[,colnames(X)%in%targets$target]
X=lapply(unique(effects$tx),function(x) X[substr(rownames(X),1,7)==x,])
names(X)=unique(effects$tx)
#duplicate effects per cell line to account for replicates
Y=lapply(unique(effects$tx),function(x) as.numeric(sapply(effects$predicted_effect[effects$tx==x],rep,2)))
names(Y)=unique(effects$tx)


###########################fit a model per tx#####################################
fit=lapply(1:32,function(x) ridge(X[[x]],Y[[x]]))

###########################get coefficients per gene per tx#######################
coefs=lapply(fit,function(x) as.matrix(coef(x$finalModel,x$bestTune$lambda)))
coefs=do.call(cbind,coefs)
colnames(coefs)=names(X)
#ignore intercept, which is 0 by definition
coefs=coefs[2:nrow(coefs),]
write.table(coefs,"coefficients_z-matrix_tx.tsv",sep='\t',quote=F)
