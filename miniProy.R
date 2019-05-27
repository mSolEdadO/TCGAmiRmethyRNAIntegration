library(biomaRt)  
library(genefu)
library(mixOmics)
load("data/porSubti.RData")

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
	filters = "ensembl_gene_id", values=rownames(concatenadas$normal[[3]]), mart=mart)

M<-pam50$centroids
PAM50genes=myannot[myannot$ensembl_gene_id%in%rownames(PAM50),]
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="ORC6",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NUF2",])
PAM50genes=rbind(PAM50genes,myannot[myannot$hgnc_symbol=="NDC80",])
write.table(PAM50genes,"data/PAM50genes.tsv",sep='\t',quote=F,row.names=F)

expr=concatenadas$normal[[3]][rownames(concatenadas$normal[[3]])%in%PAM50genes$ensembl_gene_id,]
PAM50genes=PAM50genes[order(match(PAM50genes$ensembl_gene_id,rownames(expr))),]
rownames(expr)=PAM50genes$hgnc_symbol

which(rownames(expr)=="ESR1")
#9
train = sample (1:ncol(expr) ncol(x)/2)
test =(-train)
predictors=t(do.call(rbind,concatenadas$normal[1:2]))
grid =10^ seq (10 , -2 , length =100)
#λ from 1010 to 10−2 covers the full scenarios from the null model containing only the intercept, to the least squares fit
lasso.ESR1=glmnet(x=predictors[train,],
#glmnet() function standardizes the variables so that they are on the same scale
	y=expr[9,train],
	alpha=0.5,
	lambda=grid)#larger λ: smaller coefficient 

cv.ESR1=cv.glmnet(x=predictors[train,],y=expr[9,train],alpha=0.5,nfolds=ncol(expr))
plot(cv.ESR1)
mean((predict(lasso.ESR1,s=0,newx=predictors[test,])-expr[9,test]))
#[1] -0.2096071
predi.ESR1=predict(lasso.ESR1,s=cv.ESR1$lambda.min,newx=predictors[test,])
mean((predi.ESR1-expr[9,test]))
#[1] -0.06095145
out = glmnet (x=predictors,y=expr[9,],alpha=0.5)
coef.out = predict(out, type="coefficients", s=cv.ESR1$lambda.min)

temp=coef(out)[,which(colSums(as.matrix(coef(out))>0)>1)]
temp=temp[unique(which(as.matrix(temp)>0,arr.ind=T)[,1]),]
