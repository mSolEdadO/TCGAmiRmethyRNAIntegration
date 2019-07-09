library(biomaRt)

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
enriquecimiento=lapply(sets,function(y) sapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_KEGG,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 10, 
 					pAdjustMethod = "fdr"))
 vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 vias=sapply(rownames(vias),function(x) substr(x,4,nchar(x)))
 vias=unlist(diccionario_kegg[names(diccionario_kegg)%in%vias])
return(vias)}))
sapply(enriquecimientoKEGG,function(x) table(sapply(x,class)))

#GO molecular function enrichment
enriquecimientoMF=lapply(sets,function(y) sapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_MF,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 10, 
 					pAdjustMethod = "fdr"))
 vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 vias=sapply(vias,function(x) Term(x))
return(vias)}))
i=lapply(enriquecimientoMF,function(x) which(sapply(x,length)>0))
enriquecimientoMF=lapply(1:5,function(x) 
	cbind(as.character(unique(modelos[[x]]$pam50))[i[[x]]],paste(enriquecimientoMF[[x]][i[[x]]],sep=',')))
enriquecimientoMF=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],enriquecimientoMF[[x]])))
enriquecimientoMF[,2]=sapply(enriquecimientoMF[,2],function(x) myannot$hgnc_symbol[myannot$ensembl_gene_id==x])

#GO celular component enrichment
enriquecimientoCC=lapply(sets,function(y) sapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_CC,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 10, 
 					pAdjustMethod = "fdr"))
 vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 vias=sapply(vias,function(x) Term(x))
return(vias)}))
i=lapply(enriquecimientoCC,function(x) which(sapply(x,length)>0))
enriquecimientoCC=lapply(1:5,function(x) 
	cbind(as.character(unique(modelos[[x]]$pam50))[i[[x]]],paste(enriquecimientoCC[[x]][i[[x]]],sep=',')))
enriquecimientoCC=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],enriquecimientoCC[[x]])))
enriquecimientoCC[,2]=sapply(enriquecimientoCC[,2],function(x) myannot$hgnc_symbol[myannot$ensembl_gene_id==x])

#GO bio process enrichment
enriquecimientoBP=lapply(sets,function(y) sapply(y,function(x) {
 vias=as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_BP,
 					universe= universo, 
 					hits=as.character(x),
 					minGeneSetSize = 10, 
 					pAdjustMethod = "fdr"))
 vias=rownames(vias)[vias$Adjusted.Pvalue < 0.01]
 vias=sapply(vias,function(x) Term(x))
return(vias)}))
i=lapply(enriquecimientoBP,function(x) which(sapply(x,length)>0))
enriquecimientoBP=lapply(1:5,function(x) 
	cbind(as.character(unique(modelos[[x]]$pam50))[i[[x]]],paste(enriquecimientoBP[[x]][i[[x]]],sep=',')))
enriquecimientoBP=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],enriquecimientoBP[[x]])))
enriquecimientoBP[,2]=sapply(enriquecimientoBP[,2],function(x) myannot$hgnc_symbol[myannot$ensembl_gene_id==x])

enriched=rbind(enriquecimientoCC,enriquecimientoMF,enriquecimientoBP)
colnames(enriched)=c("Subtype","gen","GO")
enriched=enriched[order(enriched[,2]),]
write.table(enriched,"enrichment.tsv",sep='\t',quote=F,row.names=F)
