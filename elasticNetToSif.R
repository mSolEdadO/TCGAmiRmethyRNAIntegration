library(caret)
library(gplots)
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

modelsToSif=function(subti){
f=files[grep(subti,files)]
	modelos=sapply(1:length(f),function(x) {
	gen=lapply(strsplit(gsub("parallel-caret/","",f[x]),".",fixed=T),function(y) y[1])
	print(paste("Loading",gen,sep=" "))	
	modelo=loadRData(f[x])
	predictores=predictors(x)
	transcri=predictores[grep("ENSG",predictores)]
	mirs=predictores[grep("hsa",predictores)]
	cpgs=predictores[grep("hsa|ENSG",x,perl=T,invert=T)]
    print("Getting quality table")
	qualy=modelo$results
	qualy=cbind(gen,qualy,length(cpgs),length(transcri),length(mirs))
	qualy=qualy[qualy$alpha==modelo$bestTune$alpha&qualy$lambda==modelo$bestTune$lambda,]
	print("Building sif")
	coefis=as.array(coef(modelo$finalModel,modelo$bestTune$lambda))
	coefis=cbind(gen,rownames(coefis),coefis)
	coefis[as.numeric(coefis[,3])!=0,]
return(list(quality=qualy,sif=coefis))})}

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
f=gsub("parallel-caret/eigen0.5/ENSG00000011426.LumA.RData","",gsub(".LumA.RData","",file))
predictores=predictors(modelo)
qualy=modelo$results
qualy=cbind(f,qualy,length(grep("hsa|ENSG",predictores,perl=T,invert=T)),length(grep("ENSG",predictores)),length(grep("hsa",predictores)))
qualy=qualy[qualy$alpha==modelo$bestTune$alpha&qualy$lambda==modelo$bestTune$lambda,]
sif=as.array(coef(modelo$finalModel,modelo$bestTune$lambda))
rm(modelo);gc()
sif=cbind(f,rownames(sif),sif)
sif=sif[sif[,3]>0,]
return(list(quality=qualy,g=sif))}
