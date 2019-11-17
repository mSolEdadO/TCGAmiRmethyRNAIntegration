#load needed libraries
library(data.table)
library(caret)
#load & pimp data
subtype=fread("Her2")
names=subtype$V1
subtype$V1=NULL
i=round(ncol(subtype)*0.8)
subtype=t(as.matrix(subtype))
colnames(subtype)=names
training=subtype[1:i,]
testing=subtype[(i+1):nrow(subtype),]

#define used function
rmse_omic=function(gen,omic){
 model=interacs[interacs$pam50==gen,]
 subset=testing[,colnames(testing)%in%model$predictor]
 if(nrow(model)==1){
 	return(RMSE(colSums(t(subset)*model$coef),testing[,colnames(testing)==gen]))}
 subset=subset[,order(match(colnames(subset),model$predictor))]
 i=which(substr(colnames(subset),1,1)==omic)
 if(omic=="all"){i=seq(1,ncol(subset))}
 if(length(i)==0){return(NA)}
 predis=colSums(t(subset[,i])*model$coef[i])
 obs=testing[,colnames(testing)==gen]
 RMSE(predis,obs)}

#use function
omicsContri=t(pbsapply(pam50$ensembl_gene_id,function(y) 
              sapply(c("c","E","h","all"),function(x) rmse_omic("Her2",y,x))))
rownames(omicsContri)=pam50$ensembl_gene_id
boxplot(omicsContri,y="log")
p.adjust(sapply(1:3,function(x)
                    wilcox.test(omicsContri[,4],
                                omicsContri[,x],
                                paired=T,
                                alternative="less")$p.val),
"fdr")
#[1] 1.358028e-08 2.682722e-05 8.492150e-08
table(apply(omicsContri,1,which.min))
 1  2  3  4 
 7  8  5 30 


