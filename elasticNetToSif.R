library(caret)
library(gplots)
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

modelsToSif=function(subti){
f=files[grep(subti,files)]
	modelos=lapply(1:length(f),function(x) {
	gen=sapply(strsplit(gsub("parallel-caret/","",f[x]),".",fixed=T),function(y) y[1])
	print(paste("Loading",gen,sep=" "))	
	load(f[x])
	predictores=predictors(model)
	transcri=predictores[grep("ENSG",predictores)]
	mirs=predictores[grep("hsa",predictores)]
	cpgs=predictores[grep("hsa|ENSG",predictores,perl=T,invert=T)]
    print("Getting quality table")
	qualy=model$results
	qualy=cbind(gen,qualy,length(cpgs),length(transcri),length(mirs))
	qualy=qualy[qualy$alpha==model$bestTune$alpha&qualy$lambda==model$bestTune$lambda,]
	print("Building sif")
	coefis=as.array(coef(model$finalModel,model$bestTune$lambda))
	coefis=cbind(gen,rownames(coefis),coefis)
	coefis=coefis[as.numeric(coefis[,3])!=0,]
return(list(quality=qualy,sif=coefis))})}

#list model files
files=list.files("parallel-caret/eigenscaled/",full.names=T)
files=files[grep("RD",files)]
files=gsub("//","/",files)

listos=c("Basal","Her2","LumB","normal")#subtype models I want
resultados=lapply(listos,modelsToSif)
names(resultados)=listos
sifs=lapply(resultados,function(x) x$sif)
quality=lapply(resultados,function(x) x$quality)
save(sifs,quality,file="modelos.RData")
#####################################
#how many predictors are selected por pam50 gene per subtype	       
sifs=lapply(sifs,function(x) x[x[,2]!="(Intercept)",])
totales=do.call(rbind,lapply(1:5,function(x) cbind(names(totales)[x],totales[[x]])))
totales=cbind(totales,rownames(totales))
g=graph.edgelist(totales[,c(1,3)],directed=F)
E(g)$weight=as.numeric(totales[,2])
g=t(as.matrix(g[unique(totales[,1]),unique(totales[,3])]))
g=g[order(match(rownames(g),pam50$ensembl_gene_id)),]
pdf("predictoresTotales.pdf")
 heatmap.2(g,trace="none",col=rev(heat.colors(100)),colRow=rainbow(4)[pam50$class],Rowv=F,Colv=F,na.color="gray",main="Predictores por gen",labRow=pam50$hgnc_symbol)
dev.off()
#####################################
#patters of expression/methylation of selected predictors per subtype
i=unique(unlist(lapply(sifs,function(x) unique(x[,2]))))
i=which(rownames(concatenadas$LumA)%in%i)
selec=do.call(cbind,lapply(concatenadas,function(x) x[i,]))
pdf("selectedHeat.pdf")
heatmap.2(selec,col=rev(heat.colors(74)),scale="none",Colv=F,trace="none",ColSideColors=rainbow(5)[subtis],symm=F,symkey=F,labCol=NA,dendrogram="row",breaks=col,labRow=NA,RowSideColors=c("cornflowerblue","brown1")[rowColors])
dev.off()
