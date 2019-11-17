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
rmse_omic=function(data,gen,omic){
 model=interacs[interacs$subtype==data&interacs$pam50==gen,]
 predis=model$predictor[substr(model$predictor,1,1)==omic]
 coefs=unlist(model$coef[substr(model$predictor,1,1)==omic])
 if(omic=="all"){
 	predis=model$predictor
 	coefs=unlist(model$coef)}
 if(length(predis)==0){return(NA)}
 tofit=paste(gen,paste(predis,collapse="+"),sep="~0+")
 fit=lm(tofit,
 		data=as.data.frame(training,stringsAsFactors=F))
 fit$coefficients=coefs
RMSE(predict(fit,as.data.frame(testing[,colnames(testing)!=gen],stringsAsFactors=F)),
testing[,colnames(testing)==gen])}

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
#[1] 2.539086e-06 6.672246e-01 8.142789e-02
table(apply(omicsContri,1,which.min))


