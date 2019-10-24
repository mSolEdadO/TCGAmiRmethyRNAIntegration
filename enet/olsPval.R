library(data.table)
library(pbapply)
#needed function
pvals=function(pam50){
  fit=lm(formula=paste(pam50,paste(interacs$predictor[interacs$pam50==pam50],collapse="+"),sep="~0+"),
  			data=training)
  stats=summary(fit)
  stats=list(predictors.p=stats$coefficients[,c(1,4)],#pval & coeff per predictor
#			 r.squared=stats$r.squared,
#			 adj.r.squared=stats$adj.r.squared,
#			 fstatistic=stats$fstatistic,
			 model.p=pf(stats$fstatistic[1],stats$fstatistic[2],stats$fstatistic[3],lower.tail=F))#model pval
return(stats)}
#load data
subtipo="Her2"
subtipo=fread(paste("../",subtipo,sep=""),sep='\t')
nombres=subtipo$V1
subtipo$V1=NULL
i=round(ncol(subtipo)*0.8)
training=as.data.frame(t(subtipo)[,1:i])#training subset should not be used for testing
colnames(training)=nombres
interacs=fread("slctdPrdctrs0.tsv")
interacs=interacs[interacs$subtype=="Her2"&interacs$predictor!=1,]

#significance values for Her2-PAM50 predictors
signi=pbsapply(unique(interacs$pam50),pvals)


