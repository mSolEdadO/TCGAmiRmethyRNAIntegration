library(ggplot2)
library(biomaRt)

#plot model quality accross subtypes
bestModels=read.table("bestModels.tsv",sep='\t')
png("modelsDevratio.png")
 ggplot(bestModels,aes(x=dev.ratio))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+ggtitle("Dev.ratio")
dev.off()

#build sifs for non null models 
bestModels=bestModels[bestModels$dev.ratio>0,]
files=apply(bestModels,1,function(x) paste(x[2],x[1],"coefs",sep='.'))
coefs=sapply(files,readLines)
names(coefs)=gsub(".coefs","",files)
coefs=sapply(coefs,function(x) t(do.call(cbind,strsplit(x,"\t"))))
coefs=do.call(rbind,lapply(1:length(files),function(x) cbind(names(coefs)[x],coefs[[x]])))
coefs=coefs[coefs[,2]!="x",]
coefs=coefs[coefs[,2]!="(Intercept)",]
coefs=coefs[order(coefs[,1]),]
#hgnc symbols instead of ensembl_ids
temp1=table(sapply(strsplit(coefs[,1],".",fixed=T),function(x) x[1]))
#los 50 tienen modelo en algun subtipo
length(temp1)
[1] 50
ids=unlist(sapply(1:49,function(x) 
	rep(as.character(pam50$hgnc_symbol)[pam50$ensembl_gene_id==names(temp1)[x]],temp1[x])))
temp=cbind(ids,coefs)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://apr2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
filters = "ensembl_gene_id", values=unique(as.character(temp[,3])),mart=mart)
temp=temp[order(temp[,3]),]
temp1=table(temp[grep("ENSG",temp[,3]),3])
ids=unlist(sapply(1:length(temp1),function(x) 
	rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(temp1)[x]],temp1[x])))
temp[grep("ENSG",temp[,3]),3]=ids
coefs1=lapply(unique(bestModels$subtype),function(x) temp[grep(x,temp[,2]),c(1,3,4)])

sum(table(unlist(lapply(coefs1,function(x) unique(x[,1]))))==5)
#32 genes tienen modelo en todos los subtipos
#[1] 32
#las 3 omicas están en los modelos de todos los subtipos
sapply(coefs1,function(x) table(substr(x[,2],1,1)))
  Basal Her2 LumA LumB non-tumor
c  2374 2192 2112 1144      1761
E   374   73 1162  346       465
h    15    4   20   30        40
coefs1=do.call(rbind,lapply(1:5,function(x) cbind(names(coefs1)[x],coefs1[[x]])))
colnames(coefs1)=c("subtype","pam50","predictor","coef")
write.table(coefs1,"slctdPrdctrs.tsv",sep='\t',quote=F,row.names=F)