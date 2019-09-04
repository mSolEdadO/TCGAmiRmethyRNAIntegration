library(ggplot2)
library(biomaRt)

#plot model quality accross subtypes
bestModels=read.table("bestModels.tsv",sep='\t')
png("modelsDevratio.png")
 ggplot(bestModels,aes(x=dev.ratio))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+ggtitle("Dev.ratio")
dev.off()

#build sifs for non null models 
files=apply(bestModels,1,function(x) paste(x[2],x[1],"coefs",sep='.'))
coefs=sapply(files,readLines)
names(coefs)=gsub(".coefs","",files)
coefs=sapply(coefs,function(x) t(do.call(cbind,strsplit(x,"\t"))))
coefs=do.call(rbind,lapply(1:195,function(x) cbind(names(coefs)[x],coefs[[x]])))
coefs=coefs[coefs[,2]!="x",]
coefs=coefs[coefs[,2]!="(Intercept)",]
coefs=coefs[order(coefs[,1]),]
#hgnc symbols instead of ensembl_ids
temp1=table(sapply(strsplit(coefs[,1],".",fixed=T),function(x) x[1]))
#49 de los 50 tienen modelo en algun subtipo
length(temp1)
[1] 49
ids=unlist(sapply(1:49,function(x) 
	rep(as.character(pam50$hgnc_symbol)[pam50$ensembl_gene_id==names(temp1)[x]],temp1[x])))
temp=cbind(ids,coefs)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
filters = "ensembl_gene_id", values=unique(as.character(temp[,3])),mart=mart)
temp=temp[order(temp[,3]),]
temp1=table(temp[grep("ENSG",temp[,3]),3])
ids=unlist(sapply(1:length(temp1),function(x) 
	rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(temp1)[x]],temp1[x])))
temp[grep("ENSG",temp[,3]),3]=ids
coefs=lapply(unique(bestModels$subtype),function(x) temp[grep(x,temp[,2]),c(1,3,4)])

#no hay modelo para los 50 genes en ningún subtipo
sapply(coefs,function(x) length(unique(x[,1])))
 Basal   Her2   LumA   LumB normal 
    38     46     33     40     38 
modelados=lapply(coefs,function(x) unique(sapply(strsplit(x[,1],".",fixed=T),function(y) y[1])))
#solo 17 tienen modelo en todos
length(intersect(intersect(intersect(intersect(modelados[[1]],modelados[[2]]),modelados[[3]]),modelados[[4]]),modelados[[5]]))
[1] 17
#las 3 omicas están en los modelos de todos los subtipos
sapply(coefs,function(x) table(substr(x[,2],1,1)))
  Basal Her2  LumA  LumB normal
c  8103 8657 50589 22071  18746
E   897  480   905   647    733
h    18  107    55   147     65
lapply(1:5,function(x) 
	write.table(coefs[[x]],paste(names(coefs)[x],"sif",sep='.'),sep='\t',quote=F,row.names=F,col.names=F))
