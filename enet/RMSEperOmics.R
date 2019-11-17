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
temp=as.numeric(omicsContri)
temp=as.data.frame(cbind(c(rep("CpGs",50),rep("transcript",50),rep("miRNAs",50),rep("all",50)),temp))
colnames(temp)=c("predictor","RMSE")
ggplot(temp,aes(y=as.numeric(as.character(RMSE)),color=predictor))+
         geom_boxplot()+
         coord_trans(y="log")+
         ylab("RMSE")+
         scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))

