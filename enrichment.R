library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(HTSanalyzeR)
library("stringr")
library(GO.db)

#load models
#basal=read.table("Basal.tsv",header=T,sep='\t')
luma=read.table("LumA.tsv",header=T,sep='\t')
lumb=read.table("LumB.tsv",header=T,sep='\t')
her2=read.table("Her2.tsv",header=T,sep='\t')
normal=read.table("normal.tsv",header=T,sep='\t')
modelos=list(basal,her2,luma,lumb,normal)
names(modelos)=c("Basal","Her2","LumA","LumB","normal")

#background data
GS_KEGG<-KeggGeneSets(species = "Hs")
GS_GO_CC<-GOGeneSets(species="Hs",ontologies=c("CC"))
GS_GO_MF<-GOGeneSets(species="Hs",ontologies=c("MF"))
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    
diccionario_kegg<-as.list(KEGGPATHID2NAME)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://apr2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","mirbase_id","external_gene_name","entrezgene"), mart=mart)
universo=as.character(myannot$entrezgene[!is.na(myannot$entrezgene)])

#entrez from selected predictors
sets=lapply(modelos,function(y) lapply(unique(y$pam50),function(x)  
	myannot$entrezgene[myannot$ensembl_gene_id%in%y$predictor[y$pam50==x]|myannot$mirbase_id%in%y$predictor[y$pam50==x]]))
sets=lapply(sets,function(x) sapply(x,function(y) y[!is.na(y)]))

#Kegg enrichment
enriquecimientoKEGG=lapply(sets,function(y) lapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_KEGG,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 15, 
 					pAdjustMethod = "fdr"))
 #vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 nombres=sapply(rownames(vias),function(x) substr(x,4,nchar(x)))
 rownames(vias)=unlist(diccionario_kegg[names(diccionario_kegg)%in%nombres])
return(vias)}))
temp=lapply(1:5,function(x) do.call(rbind,lapply(1:length(nombres[[x]]),function(y) 
	cbind(nombres[[x]][y],enriquecimientoKEGG[[x]][[y]]))))
enriquecimientoKEGG=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],temp[[x]])))
colnames(enriquecimientoKEGG)[1:2]=c("subtype","gene")

#GO molecular function enrichment
enriquecimientoMF=lapply(sets,function(y) lapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_MF,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 15, 
 					pAdjustMethod = "fdr"))
 #vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 rownames(vias)=sapply(rownames(vias),Term)
return(vias)}))
temp=lapply(1:5,function(x) do.call(rbind,lapply(1:length(nombres[[x]]),function(y) 
	cbind(nombres[[x]][y],enriquecimientoMF[[x]][[y]]))))
enriquecimientoMF=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],temp[[x]])))
colnames(enriquecimientoMF)[1:2]=c("subtype","gene")

#GO celular component enrichment
enriquecimientoCC=lapply(sets,function(y) lapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_CC,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 15, 
 					pAdjustMethod = "fdr"))
 #vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 rownames(vias)=sapply(rownames(vias),Term)
return(vias)}))
temp=lapply(1:5,function(x) do.call(rbind,lapply(1:length(nombres[[x]]),function(y) 
	cbind(nombres[[x]][y],enriquecimientoCC[[x]][[y]]))))
enriquecimientoCC=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],temp[[x]])))
colnames(enriquecimientoCC)[1:2]=c("subtype","gene")

#GO bio process enrichment
enriquecimientoBP=lapply(sets,function(y) lapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_BP,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 15, 
 					pAdjustMethod = "fdr"))
 #vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 rownames(vias)=sapply(rownames(vias),Term)
return(vias)}))
temp=lapply(1:5,function(x) do.call(rbind,lapply(1:length(nombres[[x]]),function(y) 
	cbind(nombres[[x]][y],enriquecimientoBP[[x]][[y]]))))
enriquecimientoBP=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],temp[[x]])))
colnames(enriquecimientoBP)[1:2]=c("subtype","gene")


enrichm=list(KEGG=enriquecimientoKEGG,MF=enriquecimientoMF,CC=enriquecimientoCC,BP=enriquecimientoBP)
enrichm=lapply(enrichm,function(x) x[order(x$Adjusted.Pvalue,decreasing=F),])
write.table(do.call(rbind,lapply(1:4,function(x)
 cbind(names(enrichm)[x],enrichm[[x]][enrichm[[x]]$Adjusted.Pvalue<0.05,]))),"temp",quote=F,sep='\t')
