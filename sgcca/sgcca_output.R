library(tidyverse)

files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)

sets=lapply(sets,function(x) 
	filter(x,omic=="transcripts")%>%
	select(c("variable","component")))

#######################################CHECK ENRICHMENTS
library(biomaRt)
library(org.Hs.eg.db)
library(KEGGREST)
library(HTSanalyzeR)
library("stringr")#do I need this???
library(GO.db)
library(GSEABase)

#get dbs 
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	host="http://feb2021.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"entrezgene_accession"), mart=mart)
universo=as.character(myannot$entrezgene_id[!is.na(myannot$entrezgene_id)])
#[1] 29067 ids
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))#list of entrez_ids                   
#[1] 12319 BPs
diccionario_kegg <- keggList("pathway", "hsa")
#[1] 345 pathways
names(diccionario_kegg)=gsub("path:","",names(diccionario_kegg))
GS_KEGG<-sapply(names(diccionario_kegg),function(x) 
	keggGet(x)[[1]]$GENE)#list like:
#$hsa00010		#path id
#  [1] "3101"   #entrez_id                                                                                      
#  [2] "HK3; hexokinase 3 [KO:K00844] [EC:2.7.1.1]" #entrezgene_accession;...                                                  

#separate symbols from IDs
#GS_KEGG_symbol=sapply(GS_KEGG,function(x) x[grep(" ",x)])
#GS_KEGG_symbol=lapply(GS_KEGG_symbol,function(x) 
	#keep only the 1st part of data
#	sapply(strsplit(as.character(x),";"),function(y) y[1]))
GS_KEGG_id=sapply(GS_KEGG,function(x) x[grep(" ",x,invert=T)])
#17 paths have entrez_ids not in biomaRt
#319 BPs have entrez_ids not in biomaRt

#translate ensembl_ids
colnames(myannot)[1]="variable"
sets=lapply(sets,function(x) 
	merge(x,myannot,by="variable"))#this merge shows ensembl ids mapping to several entrez
#$Basal
# ENSG00000143384        574406 ADAMTSL4-AS1        
# ENSG00000143384          4170 MCL1                
#$LumA
# ENSG00000159216     100506403 LOC100506403        
# ENSG00000159216           861 RUNX1

#make diccionario... mergable
diccionario_kegg=gsub(" - Homo sapiens (human)","",
	diccionario_kegg,fixed=T)
diccionario_kegg=as.data.frame(diccionario_kegg)
diccionario_kegg$hsa=rownames(diccionario_kegg)

#########KEGG enrichment
enrichKEGG=lapply(sets,function(y)
	do.call(rbind,lapply(unique(y$component),function(x) {
	#lapply per comp coz I couldn't make it with group_map	
	resus=as.data.frame(multiHyperGeoTest(
		collectionOfGeneSets=GS_KEGG_id,
		universe=universo,
		hits=as.character(y$entrezgene_id[y$component==x]),
		pAdjustMethod="fdr"))
	resus$component=x
	resus$hsa=rownames(resus)
	resus=merge(resus,diccionario_kegg,by="hsa")
})))

#enrichment of biological processes
enrichBP=lapply(sets,function(y)
	do.call(rbind,lapply(unique(y$component),function(x) {
	#lapply per comp coz I couldn't make it with group_map	
	resus=as.data.frame(multiHyperGeoTest(
		collectionOfGeneSets=GS_GO_BP,
		universe=universo,
		hits=as.character(y$entrezgene_id[y$component==x]),
		pAdjustMethod="fdr"))
	resus$component=x
	resus$BP=Term(rownames(resus))
})))

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

enrichKEGG=lapply(enrichKEGG,function(x) x[x$Adjusted.Pvalue<0.05,])
enrichKEGG=do.call(rbind,lapply(which(sapply(enrichKEGG,nrow)>0),
	function(x) cbind(names(enrichKEGG)[x],rownames(enrichKEGG[[x]]),
		enrichKEGG[[x]])))
colnames(enrichKEGG)[1:2]=c("PC","hsa")
enrichKEGG$name=unlist(sapply(gsub("hsa","",enrichKEGG$hsa),function(x) diccionario_kegg[names(diccionario_kegg)%in%x]))
