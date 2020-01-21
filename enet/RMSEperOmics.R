#load needed libraries
library(data.table)
library(caret)
library(gplots)
library(ggsignif)
library(gridExtra)
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

#repeat per subtype!!!!
omicsContri=t(pbsapply(pam50$ensembl_gene_id,function(y) 
              sapply(c("c","E","h","all"),function(x) rmse_omic(y,x))))
rownames(omicsContri)=pam50$ensembl_gene_id
#table(apply(omicsContri,1,which.min))
# 1  2  3  4 
# 7  8  5 30 

#plot per subtype
omicsContri=list(Basal,Her2,LumB,LumA,normal)
omicsContri=lapply(omicsContri,function(x) 
 cbind(unlist(x),as.character(sapply(c("CpG","gene","miRNA","all"),rep,50))))
plots=lapply(c(2,3,1,4,5),function(x) 
             ggplot(omicsContri[[x]],aes(x=as.numeric(as.character(RMSE))))+
             geom_density(aes(group=predictor,color=predictor,fill=predictor,y=..scaled..),alpha=0.3)+
             ylab("scaled frequency")+xlab("testing RMSE")+
             scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+
             scale_fill_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+xlim(1,1e+5)+
             scale_x_continuous(trans='log10',breaks=c(100,1000,10000))+ggtitle(names(omicsContri)[x]))
png("omicsContrib.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],ncol=2)
dev.off()
lapply(subtypes,function(y) ks.test(as.numeric(as.character(y$x[y$predictor=="CpG"])),as.numeric(as.character(y$x[y$predictor=="miRNA"]))))
#Warning messages:
#1: In ks.test(as.numeric(as.character(y$x[y$predictor == "CpG"])),  :
#  cannot compute exact p-value with ties

#plot per omic
omics=lapply(unique(omicsContri[[1]][,2]),function(y) 
 as.data.frame(cbind(as.numeric(sapply(omicsContri,function(x) x[x[,2]==y,1])),
                     as.character(sapply(c("Basal","Her2","LumB","LumA","normal"),rep,50)))))
colnames(omics[[1]])=c("RMSE","subtype")
plots=lapply(1:4,function(x) 
            ggplot(omics[[x]],aes(x=as.numeric(as.character(RMSE))))+
            geom_density(aes(group=subtype,color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
            ylab("scaled frequency")+xlab("testing RMSE")+
            scale_x_continuous(trans='log10',breaks=c(100,1000,10000))+ggtitle(names(omics)[x]))
png("SubtypeOmicsContrib.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],ncol=2)
dev.off()
