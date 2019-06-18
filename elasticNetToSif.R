library(caret)
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

#load files
files=list.files()
files=files[grep("RD",files)]
files=gsub("//","/",files)
temp=files[grep("LumA",files)]
modelos=list()
modelos=lapply(1:length(herFiles),function(x) modelos[[x]]=loadRData(temp[x]))
temp=gsub(".LumA.RData","",gsub("parallel-caret/eigenscaled/","",temp))

predictores=lapply(modelos,function(x) predictors(x))
transcri=sapply(predictores,function(x) x[grep("ENSG",x)])
mirs=lapply(predictores,function(x) x[grep("hsa",x)])
cpgs=lapply(predictores,function(x) x[grep("hsa|ENSG",x,perl=T,invert=T)])

qualy=lapply(modelos,function(x) x$results)
qualy=lapply(1:45,function(x) qualy[[x]][qualy[[x]]$alpha==modelos[[x]]$bestTune$alpha&qualy[[x]]$lambda==modelos[[x]]$bestTune$lambda,])
qualy=cbind(temp,do.call(rbind,qualy))
lumaQualy=cbind(qualy,sapply(cpgs,length),sapply(transcri,length),sapply(mirs,length))

temp=lumaQualy[,c(1,10:12)]
lala=lumbQualy[!lumbQualy$ensembl%in%temp$ensembl,c(1,10:12)]
lala$CpGs=NA
lala$transcri=NA
lala$miR=NA
temp=rbind(temp,lala)
temp=temp[order(match(temp$ensembl,lumbQualy$ensembl)),]

coefis=lapply(1:45,function(x) as.array(coef(modelos[[x]]$finalModel,modelos[[x]]$bestTune$lambda)))
coefis=lapply(coefis,function(x) x[x>0,])
coefis=lapply(1:45,function(x) cbind(temp[x],coefis[[x]]))
lumaSif=do.call(rbind,coefis)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
mannot=getBM(attributes = c("chromosome_name","ensembl_gene_id","hgnc_symbol"),filters = "ensembl_gen
e_id", values=normalQualy$ensembl, mart=mart)
diferencias=list(CpGs=cbind(unlist(normalQualy$CpGs),unlist(lumbQualy$CpGs),unlist(basalQualy$CpGs),unlist(h
erQualy[,10])),transcri=cbind(unlist(normalQualy$transcri),unlist(lumbQualy$transcri),unlist(basalQua
ly$transcri),unlist(herQualy[,11])),mir=cbind(unlist(normalQualy$miR),unlist(lumbQualy$miR),unlist(ba
salQualy$miR),unlist(herQualy[,12])))
rownames(diferencias$CpGs)=mannot$hgnc_symbol

totales=diferencias$CpGs+diferencias$transcri+diferencias$mir
proporciones=lapply(diferencias,function(x) x/totales)
proporciones[[2]][is.na(proporciones[[2]])]=0
proporciones[[2]][is.na(totales)]=NA

pdf("proporPredictores.pdf")
 heatmap.2(proporciones[[1]],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[temp$class],Rowv=F,Colv=F,na.color="gray",main="CpGs/total")
 heatmap.2(proporciones[[2]],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[temp$class],Rowv=F,Colv=F,na.color="gray",main="transcri/total")
 heatmap.2(proporciones[[3]],dendrogram="none",trace="none",col=c("white",rev(heat.colors(100))),colRow=rainbow(4)[temp$class],Rowv=F,Colv=F,na.color="gray",main="miR/total")
dev.off()

sifs=list(normal=normalSif,basal=basalSif,her2=herSif,lumA=lumaSif,lumB=lumBsif)
quality=list(normal=normalQualy,basal=basalQualy,her2=herQualy,lumA=lumaQualy,lumB=lumbQualy)
save(sifs,quality,file="modelosEigenScale.RData")

nodos=lapply(sifs,function(x) unique(c(rownames(x),x[,1])))
nodos=lapply(nodos,function(x) x[x!="(Intercept)"])

