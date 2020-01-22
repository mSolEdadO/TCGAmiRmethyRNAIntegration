###########################libraries###############################################
library(data.table)
library(caret)
library(doParallel)
#registerDoParallel(cores=2)
#ugly function to read perturbation matrices
load_pertur=function(file){
 temp=fread(file)
 names=temp$hgnc_symbol
 temp=as.matrix(temp[,2:ncol(temp)])
 rownames(temp)=names
return(temp)}
#ridge regression function
ridge=function(expression,inhibition){
	#parameters to fit: lambda=penalization chosen by cv, alpha=0 since this is a ridge regression
	coefGrid <-  expand.grid(lambda=seq (0,50, length =100),#length is usually 100
			alpha=0)
	set.seed(123)
	#cv parameters
	trainCtrl <- trainControl("repeatedcv",
			 number = 3, #min number of folds, 3 since n is small 
			 repeats=10,#the more the better li
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
files=list.files("/labs/csbig/DREAMCHALLENGE_CTD2/raw/",full.names=T)
#load response data per cell line
response=lapply(files[grep("Response",files)],fread)
names(response)=sapply(strsplit(sapply(strsplit(files[grep("Response",files)],"/"),
				function(x) x[7]),"-"),function(y) y[1])
#load perturbation data per cell line 
perturbation=lapply(files[grep("RNA",files)],load_pertur)
perturbation=do.call(cbind,perturbation)

###########################formating data################################
#transform the list of perturbations per cell line into a list of perturbation per tx
tx=unique(sapply(strsplit(colnames(perturbation),"_"),function(x) paste(x[1:2],collapse='_')))
perturbation=lapply(tx,function(x) perturbation[,grep(x,colnames(perturbation))])
names(perturbation)=tx
concentrations=read.table("/labs/csbig/DREAMCHALLENGE_CTD2/results/experimentAnnotation_df.tsv",header=T)
#get the nearest response to the concetrations tested, this is the Y to be modeled
Y=lapply(levels(concentrations$compound_id),function(x) sapply(levels(concentrations$cell_line),function(y) response[[y]][[x]][!is.na(response[[y]][[x]])][which.min(unique(concentrations$concentration[concentrations$compound_id==x&concentrations$cell_line==y])-10^response[[y]]$dose_log10_uM[!is.na(response[[y]][[x]])])]))
names(Y)=levels(concentrations$compound_id)
#duplicate to account for replicates
Y=lapply(Y,function(x) c(x,x))
Y=lapply(Y,function(x) x[order(names(x))])
#perturbation has the X to model Y
perturbation=perturbation[order(match(names(perturbation),names(Y)))]
#features should be in columns
perturbation=lapply(perturbation,t)
#get the cell lines
cell.line=lapply(perturbation,function(x) 
	sapply(strsplit(rownames(x),"_"),function(y) y[5]))
#X to model Y taking cell lines as another feature
tomodel=lapply(1:32,function(x) 
	model.matrix(~0+apply(perturbation[[x]],2,function(x) x)+factor(cell.line[[x]])))

###########################fit a model per tx###########################						      
fit=lapply(1:32,function(x) ridge(tomodel[[x]],Y[[x]]))
#coefs=as.matrix(coef(model$finalModel, model$bestTune$lambda))
