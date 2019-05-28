library(biomaRt)  
library(genefu)
library(caret)
load("data/porSubti.RData")
load("data/MFAmultio.RData")

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
	filters = "ensembl_gene_id", values=rownames(concatenadas$normal[[3]]), mart=mart)

M<-pam50$centroids
PAM50genes=myannot[myannot$ensembl_gene_id%in%rownames(PAM50),]
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="ORC6",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NUF2",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NDC80",])
write.table(PAM50genes,"data/PAM50genes.tsv",sep='\t',quote=F,row.names=F)

normalEigenScaled=MFAmultio$normal$global.pca$call$X
#preProcess
#near-zero-variance” predictors may need to be identified and eliminated prior to modeling
nzv <- nearZeroVar(normalEigenScaled, saveMetrics= TRUE)
table(nzv$zeroVar)
# FALSE 
#398925 !!!!!!!!!!!!!!!!!!!
#i=which(colnames(normalEigenScaled)==PAM50genes$ensembl_gene_id[PAM50genes$hgnc_symbol=="ESR1"])
i=which(colnames(normalEigenScaled)==pam50$ensembl_gene_id[pam50$hgnc_symbol=="ESR1"])
y=normalEigenScaled[,i] 
normalEigenScaled=normalEigenScaled[,1:ncol(normalEigenScaled)!=i]

set.seed(123)
train = sample (1:nrow(normalEigenScaled),2*nrow(normalEigenScaled)/3)
test =(-train)
grid=seq(0,1,length.out=100)#should cover the full scenarios from the null model containing 
#only the intercept, to the least squares fit
cv.ESR1=cv.glmnet(x=as.matrix(normalEigenScaled[train,]),
	y=y[train],
	standardize=T,#data is alredy tranformed
	alpha=0.5,
	nfolds=5,
	lambda=grid,#larger λ: smaller coefficient
	type.gaussian = "naive")#p>500
#glmnet chooses the lambda sequence such that the number of nonzero coefficients ranges from 0 
#to p, where p is the total number of features
plot(cv.ESR1)
out = glmnet (x=as.matrix(normalEigenScaled),y=y,alpha=0.5,lambda=grid)
coef.out = predict(out, type="coefficients", s=cv.ESR1$lambda.min) #min vs 1se
#lambda.1se represents the λ that was simpler than the best model(lambda.min), but has error 
#within 1 standard error of the best model, so it cannot be distinguished from the best model 
mean((predict(lasso.ESR1,s=0,newx=predictors[test,])-expr[9,test]))
predi.ESR1=predict(lasso.ESR1,s=cv.ESR1$lambda.min,newx=predictors[test,])
mean((predi.ESR1-expr[9,test]))

#####################################################

set.seed(123)
model <- train(
  y~., data = normalEigScldPCA, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Best tuning parameter
model$bestTune
coef(model$finalModel, model$bestTune$lambda)