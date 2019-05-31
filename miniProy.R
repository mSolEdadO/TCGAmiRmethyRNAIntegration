library(biomaRt)  
library(genefu)
library(caret)
load("data/porSubti.RData")
load("data/MFAmultio.RData")

#get pam50 genes ids
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
	filters = "ensembl_gene_id", values=rownames(concatenadas$normal[[3]]), mart=mart)
M<-pam50$centroids#genefu pam50 is on hgnc symbols
PAM50genes=myannot[myannot$ensembl_gene_id%in%rownames(PAM50),]
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="ORC6",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NUF2",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NDC80",])
write.table(PAM50genes,"data/PAM50genes.tsv",sep='\t',quote=F,row.names=F)

#############################
#try for a single "subtype"
normalEigenScaled=MFAmultio$normal$global.pca$call$X
#preProcess
#near-zero-variance predictors may need to be identified and eliminated prior to modeling
nzv <- nearZeroVar(normalEigenScaled, saveMetrics= TRUE)
table(nzv$zeroVar)
# FALSE 
#398925 !!!!!!!!!!!!!!!!!!!
#i=which(colnames(normalEigenScaled)==PAM50genes$ensembl_gene_id[PAM50genes$hgnc_symbol=="ESR1"])
i=which(colnames(normalEigenScaled)==pam50$ensembl_gene_id[pam50$hgnc_symbol=="ESR1"])
y=normalEigenScaled[,i] 
normalEigenScaled=normalEigenScaled[,1:ncol(normalEigenScaled)!=i]

#cross validation
grid =10^ seq (10 , -2 , length =100)#cover the full scenarios from the null model containing 
#only the intercept, to the least squares fit
cv.ESR1=cv.glmnet(x=as.matrix(normalEigenScaled), #cv creates the testing subsets,
	y=y,
	standardize=F,#data is alredy tranformed
	alpha=0.5,
	nfolds=10,
	lambda=grid,#larger λ: smaller coefficient
	type.gaussian = "naive")#p>500
plot(cv.ESR1)
#model
out = glmnet (x=as.matrix(normalEigenScaled),y=y,alpha=0.5,lambda=grid)
#get coefficients
coef.out = predict(out, type="coefficients", s=cv.ESR1$lambda.min) #min vs 1se
#lambda.1se represents the λ that was simpler than the best model(lambda.min), but has error 
#within 1 standard error of the best model, so it cannot be distinguished from the best model 
########################################################
#repeat for all subtypes
library(parallel)

modelo=function(y,x){
cvmodel=cv.glmnet(x=as.matrix(x),y=y,standardize=F,alpha=0.5,nfolds=10,lambda=grid,type.gaussian = "naive")#p>500
#plot(cvmodel)
out = glmnet (x=as.matrix(x),y=y,standardize=F,alpha=0.5,lambda=grid,type.gaussian = "naive")
coef.out = predict(out, type="coefficients", s=cvmodel$lambda.min) #min vs 1se
colnames(coef.out)=cvmodel$lambda.min
return(coef.out)
}

no_cores <- detectCores() - 2  
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(glmnet)})
#clusterEvalQ(cl,{library(pbapply)})
clusterExport(cl,c("pam50","modelo","grid"))
#model per subtype
coeficientes=parLapply(cl,Eigenscaled,function(s) 
	sapply(pam50$ensembl_gene_id,function(g) 
		modelo(y=s[,colnames(s)==g],x=s[,colnames(s)!=g])))
stopCluster(cl)  
coeficientes=lapply(coeficientes,function(x) sapply(x,function(y) as.matrix(y[which(y>0),])))
coeficientes=lapply(coeficientes,function(y) do.call(rbind,sapply(1:length(y),function(x) 
	cbind(names(y)[x],rownames(y[[x]]),y[[x]]))))

#Annotation
#downloaded last genome build CpG annotation from:
#http://zwdzwd.io/InfiniumAnnotation/current/hm450/hm450.hg38.manifest.rds
#Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes
#Wanding Zhou Peter W. Laird Hui Shen
methy=read.table("../Downloads/hm450.hg38.manifest.tsv.gz",sep='\t',header=T)
nodos=unique(unlist(normal[,1:2]))
methy=methy[methy$probeID%in%nodos[grep("cg",nodos)],c(1:3,5,21)]
#mannot=getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"),
#	filters = "ensembl_gene_id", values=nodos[grep("ENS",nodos)], mart=mart)
mannot=getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol"),
	filters = "hgnc_symbol", values=nodos, mart=mart)
mirannot=getBM(attributes = c("chromosome_name","start_position","end_position","mirbase_id"),
 mart=mart,filter="mirbase_id",values=nodos)
