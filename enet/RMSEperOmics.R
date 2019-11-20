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
omicsContri=rbind(lumA,lumB,basal,her2,normal)
ks.test(omicsContri[,1],omicsContri[,2],alternative="less")
ks.test(omicsContri[,1],omicsContri[,3],alternative="less")
ks.test(omicsContri[,1],omicsContri[,4],alternative="less")
ks.test(omicsContri[,4],omicsContri[,1],alternative="greater")
ks.test(omicsContri[,4],omicsContri[,2],alternative="greater")
ks.test(omicsContri[,4],omicsContri[,3],alternative="greater")

RMSE=as.numeric(unlist(omicsContri))
predictor=c(rep("CpG",250),rep("transcript",250),rep("miRNA",250),rep("mix",250))
omicsContri=as.data.frame(cbind(RMSE,predictor),stringsAsFactors=F)
boxes=ggplot(omicsContri,aes(y=as.numeric(RMSE),color=predictor,x=predictor))+
      geom_boxplot()+ylab("testing RMSE")+ylim(0,280000)+
      scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+
      geom_signif(test="ks.test",
                  comparisons=list(c("all","CpG"),c("all","miRNA"),c("all","transcript"),c("CpG","miRNA"),c("CpG","transcript")),
                  map_signif_level=T,y_position=c(185000,215000,265000,190000,245000))
densi=ggplot(omicsContri,aes(x=as.numeric(RMSE)))+
      geom_density(aes(group=predictor,color=predictor,fill=predictor,y=..scaled..),alpha=0.3)+
      ylab("scaled frequency")+xlab("testing RMSE")+scale_x_continuous(trans='log10')+
      scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+
      scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+
png("OmicsContrib.png")
 grid.arrange(boxes,densi,nrow=2)
dev.off()

plots=lapply(c(2,3,1,4,5),function(x) ggplot(omicsContri[[x]],aes(x=as.numeric(as.character(RMSE))))+geom_density(aes(group=predictor,color=predictor,fill=predictor,y=..scaled..),alpha=0.3)+ylab("scaled frequency")+xlab("testing RMSE")+scale_color_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+scale_fill_manual(values=c("firebrick1","#999999", "#E69F00", "#56B4E9"))+xlim(1,1e+5)+scale_x_continuous(trans='log10',breaks=c(100,1000,10000))+ggtitle(names(omicsContri)[x]))
png("omicsContrib.png")
> grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],ncol=2)
dev.off()