library(data.table)
library(caret)
library(doParallel)
#registerDoParallel(cores=2)

load_pertur=function(file){
 temp=fread(file)
 names=temp$hgnc_symbol
 temp=as.matrix(temp[,2:ncol(temp)])
 rownames(temp)=names
return(temp)}


files=list.files("/labs/csbig/DREAMCHALLENGE_CTD2/raw/",full.names=T)
response=lapply(files[grep("Response",files)],fread)
names(response)=sapply(strsplit(sapply(strsplit(files[grep("Response",files)],"/"),
				function(x) x[7]),"-"),function(y) y[1])
perturbation=lapply(files[grep("RNA",files)],load_pertur)
perturbation=do.call(cbind,perturbation)
tx=unique(sapply(strsplit(colnames(perturbation),"_"),function(x) paste(x[1:2],collapse='_')))
perturbation=lapply(tx,function(x) perturbation[,grep(x,colnames(perturbation))])
names(perturbation)=tx

concentrations=read.table("experimentAnnotation_df.tsv",header=T)
Y=lapply(levels(concentrations$compound_id),function(x) sapply(levels(concentrations$cell_line),function(y) response[[y]][[x]][!is.na(response[[y]][[x]])][which.min(unique(concentrations$concentration[concentrations$compound_id==x&concentrations$cell_line==y])-10^response[[y]]$dose_log10_uM[!is.na(response[[y]][[x]])])]))
names(Y)=levels(concentrations$compound_id)
Y=lapply(Y,function(x) c(x,x))
Y=lapply(Y,function(x) x[order(names(x))])
perturbation=perturbation[order(match(names(perturbation),names(Y)))]
perturbation=lapply(perturbation,t)
cell.line=lapply(perturbation,function(x) 
	sapply(strsplit(rownames(x),"_"),function(y) y[5]))
tomodel=lapply(1:32,function(x) 
	model.matrix(~0+apply(perturbation[[x]],2,function(x) x)+factor(cell.line[[x]])))
						      
ridge=function(expression,inhibition){
	coefGrid <-  expand.grid(lambda=seq (0,50, length =100),#length is usually 100
			alpha=0)
	set.seed(123)
	trainCtrl <- trainControl("repeatedcv",
			 number = 3, #k choose this according to n
			 repeats=10,#200 its ok, the more the better li
			 verboseIter = T,#T if fit fails,
			 allowParallel=T)
	set.seed(567)
	model <- train(y = inhibition,
	       x = expression,
	       method = "glmnet",
	       trControl = trainCtrl,
	       tuneGrid = coefGrid,
	       standardize=T,
	       intercept=F)
	return(model)}
fit=lapply(1:32,function(x) ridge(tomodel[[x]],Y[[x]]))
#coefs=as.matrix(coef(model$finalModel, model$bestTune$lambda))
