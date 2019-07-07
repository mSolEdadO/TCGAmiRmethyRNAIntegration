library(caret)
library(gplots)
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

modelsToSif=function(subti){
f=files[grep(subti,files)]
modelos=list()
print("Loading files")
modelos=lapply(1:length(f),function(x) modelos[[x]]=loadRData(f[x]))
f=sapply(strsplit(gsub("parallel-caret/eigenscaled/","",f),".",fixed=T),function(x) x[1])

print("Getting quality table")
predictores=lapply(modelos,function(x) predictors(x))
transcri=sapply(predictores,function(x) x[grep("ENSG",x)])
mirs=lapply(predictores,function(x) x[grep("hsa",x)])
cpgs=lapply(predictores,function(x) x[grep("hsa|ENSG",x,perl=T,invert=T)])
qualy=lapply(modelos,function(x) x$results)
qualy=lapply(1:length(modelos),function(x) 
qualy[[x]][qualy[[x]]$alpha==modelos[[x]]$bestTune$alpha&qualy[[x]]$lambda==modelos[[x]]$bestTune$lambda,])
qualy=cbind(f,do.call(rbind,qualy))
qualy=cbind(qualy,sapply(cpgs,length),sapply(transcri,length),sapply(mirs,length))

print("Building sif")
coefis=lapply(1:length(modelos),function(x) 
as.array(coef(modelos[[x]]$finalModel,modelos[[x]]$bestTune$lambda)))
coefis=lapply(1:length(modelos),function(x) cbind(f[x],rownames(coefis[[x]]),coefis[[x]]))
#coefis=lapply(coefis,function(x) x[x>0,])
sif=do.call(rbind,coefis)
return(list(quality=qualy,sif=sif))}

#list model files
files=list.files("parallel-caret/eigenscaled/",full.names=T)
files=files[grep("RD",files)]
files=gsub("//","/",files)

listos=c("Basal","Her2","LumB","normal")#subtype models I want
resultados=lapply(listos,modelsToSif)
names(resultados)=listos
sifs=lapply(resultados,function(x) x$sif[as.numeric(x$sif[,3])!=0,])#needed!
quality=lapply(resultados,function(x) x$quality)
save(sifs,quality,file="modelos.RData")
	       
#the number of selected predictors is the same across subtypes?
prediTot=do.call(cbind,lapply(quality,function(x) rowSums(x[,10:12])))
rownames(prediTot)=quality$Her2[,1]
geneClass=read.table("resultados/pam50.tsv")
prediTot=prediTot[order(match(rownames(prediTot),geneClass$ensembl_gene_id)),]
prediCs=do.call(cbind,lapply(quality,function(x) x[,10]))
prediTs=do.call(cbind,lapply(quality,function(x) x[,11]))
prediMs=do.call(cbind,lapply(quality,function(x) x[,12]))
rownames(prediCs)=quality$Her2[,1]
rownames(prediTs)=quality$Her2[,1]
rownames(prediMs)=quality$Her2[,1]
prediCs=prediCs[order(match(rownames(prediCs),geneClass$ensembl_gene_id)),]
prediTs=prediTs[order(match(rownames(prediTs),geneClass$ensembl_gene_id)),]
prediMs=prediMs[order(match(rownames(prediMs),geneClass$ensembl_gene_id)),]   
pdf("proporPredictores05.pdf")
 heatmap.2(prediTot[c(1:44,46:50),],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[geneClass$class],Rowv=F,Colv=F,na.color="gray",main="total",labRow=geneClass$hgnc_symbol)
 heatmap.2(prediCs[c(1:44,46:50),],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[geneClass$class],Rowv=F,Colv=F,na.color="gray",main="total",labRow=geneClass$hgnc_symbol)
 heatmap.2(prediTs[c(1:44,46:50),],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[geneClass$class],Rowv=F,Colv=F,na.color="gray",main="total",labRow=geneClass$hgnc_symbol)
 heatmap.2(prediMs[c(1:44,46:50),],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[geneClass$class],Rowv=F,Colv=F,na.color="gray",main="total",labRow=geneClass$hgnc_symbol)
dev.off()
###############################
#pa' lumA			     
modelsToSif=function(file){
modelo=loadRData(file)
f=gsub("parallel-caret/","",gsub(".LumA.RData","",file))
predictores=predictors(modelo)
qualy=modelo$results
qualy=cbind(f,qualy,length(grep("hsa|ENSG",predictores,perl=T,invert=T)),length(grep("ENSG",predictores)),length(grep("hsa",predictores)))
qualy=qualy[qualy$alpha==modelo$bestTune$alpha&qualy$lambda==modelo$bestTune$lambda,]
sif=as.array(coef(modelo$finalModel,modelo$bestTune$lambda))
rm(modelo);gc()
sif=cbind(f,rownames(sif),sif)
sif=sif[sif[,3]>0,]
return(list(quality=qualy,g=sif))}

