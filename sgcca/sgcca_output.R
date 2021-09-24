#######################################CHECK ENRICHMENTS
library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(HTSanalyzeR)
library("stringr")
library(GO.db)
library(GSEABase)

#get transcripts per component
selected=apply(final$loadings$transcripts,2,function(x) 
	rownames(final$loadings$transcripts)[x!=0])
#around 22 per pc with penalty
#around 253 per pc with penalty1
#get dbs 
GS_KEGG<-KeggGeneSets(species = "Hs")
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    
diccionario_kegg<-as.list(KEGGPATHID2NAME)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol",
	"entrezgene_id"), mart=mart)
universo=as.character(myannot$entrezgene_id[!is.na(myannot$entrezgene_id)])
#transform to entrez id
sets=lapply(selected,function(y) myannot$entrezgene_id[myannot$ensembl_gene_id%in%y])
sets=lapply(sets,function(x) x[!is.na(x)])
#enrichment of biological processes
enrichBP=lapply(sets,function(y)
	as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_BP,
 	universe= universo, 
 	hits=as.character(y),
 	minGeneSetSize = 15, 
 	pAdjustMethod = "fdr")))
#make it a data.frame
enrichBP=lapply(enrichBP,function(x) x[x$Adjusted.Pvalue<0.05,])
enrichBP=do.call(rbind,lapply(which(sapply(enrichBP,nrow)>0),
	function(x) cbind(names(enrichBP)[x],rownames(enrichBP[[x]]),
		enrichBP[[x]])))
colnames(enrichBP)[1:2]=c("PC","GO")
enrichBP$name=Term(as.character(enrichBP$GO))
fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
slim <- getOBOCollection(fl)
bpGroups=goSlim(myCollection, slim, "BP")



#net per enriched component
source("function_networkAlt.R")
g=network(final,comp=list(CpGs=7,transcripts=7,miRNAs=7),
	blocks=1:3)$gR
temp=as.data.frame(get.edgelist(g))
temp$cor=E(g)$weight

enrichKEGG=lapply(sets,function(y)
	as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_KEGG,
 	universe= universo, 
 	hits=as.character(y),
 	minGeneSetSize = 15, 
 	pAdjustMethod = "fdr")))
enrichKEGG=lapply(enrichKEGG,function(x) x[x$Adjusted.Pvalue<0.05,])
enrichKEGG=do.call(rbind,lapply(which(sapply(enrichKEGG,nrow)>0),
	function(x) cbind(names(enrichKEGG)[x],rownames(enrichKEGG[[x]]),
		enrichKEGG[[x]])))
colnames(enrichKEGG)[1:2]=c("PC","hsa")
enrichKEGG$name=unlist(sapply(gsub("hsa","",enrichKEGG$hsa),function(x) diccionario_kegg[names(diccionario_kegg)%in%x]))
