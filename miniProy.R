library(biomaRt)  
library(genefu)
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
grid =10^ seq (10 , -2 , length =100)
#λ from 1010 to 10−2 covers the full scenarios from the null model containing only the intercept, to the least squares fit
lasso.ESR1=glmnet(x=t(do.call(rbind,concatenadas$normal[1:2])),
	y=expr[9,],alpha=0.5,lambda=grid)
#glmnet() function standardizes the variables so that they are on the same scale
#larger λ: smaller coefficient 
temp=coef(lasso.ESR1)[,which(colSums(as.matrix(coef(lasso.ESR1))>0)>1)]
temp=temp[unique(which(as.matrix(temp)>0,arr.ind=T)[,1]),]
