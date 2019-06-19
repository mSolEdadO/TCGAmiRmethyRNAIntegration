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
############################################################################
#validate
library(biomaRt)
library(rentrez)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name","entrezgene"), mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]

#hgnc_symbol instead of ensemblID
interacs=unique(do.call(rbind,sifs))
interacs=interacs[order(interacs[,1]),]
interacs=interacs[interacs[,2]!="",]
interacs=interacs[interacs[,2]!="(Intercept)",]
targets=table(interacs[,1])
temp=sapply(1:50,function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(targets)[x]],targets[x]))
interacs=cbind(interacs,unlist(temp))
interacs=interacs[order(interacs[,2]),]
i=grep("ENSG",interacs[,2])
interacs=cbind(interacs,interacs[,2])
notarg=table(interacs[i,4])
temp=sapply(1:length(notarg),function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(notarg)[x]],notarg[x]))
interacs[i,4]=unlist(temp)

#how many times both entities are mentioned in the same paper
comention=apply(interacs[,3:4],1,function(x)  
	rentrez::entrez_search(db = "pubmed", term = paste(x,collapse=" "))$count)


methy=read.table("ini/hm450.hg38.manifest.tsv",sep='\t',header=T,fill=T)

#am I linking genes with their known regulating CpGs? NOP
i=grep("ENSG|hsa",interacs[,2],perl=T,invert=T)
length(i)
#[1] 3997
knownReg=sapply(i,function(x) 
	grep(interacs[x,3],methy$gene_HGNC[methy$probeID==interacs[x,4]]))
sum(sapply(knownReg,length)==0)
#[1] 3997
#are they in the same chr?
sameChr=apply(interacs[i,1:2],1,function(x) 
	sum(myannot$chromosome_name[myannot$ensembl_gene_id==x[1]]==methy$CpG_chrm[methy$probeID==x[2]]))
sum(sameChr!=0) 
#[1] 156 sólo estos están en el mismo cromosoma y están lejos ↓
apply(methy[methy$probeID%in%temp[,2],2:3]-myannot$start_position[myannot$ensembl_gene_id=="ENSG00000164611"],2,function(x) min(abs(x)))
#CpG_beg CpG_end 
#4076136 4076138 
interacs=cbind(interacs,NA)
interacs[i,5]="n"
interacs[i[sameChr!=0],5]="y"
temp=interacs[interacs[,5]=="y",1:2]#todos afectan PTTG1  
lapply(sifs,function(x) sum(x[,1]=="ENSG00000164611"&x[,2]%in%temp[,2]))
$Basal
[1] 54
$Her2
[1] 72
$LumB
[1] 30
$normal
[1] 0
