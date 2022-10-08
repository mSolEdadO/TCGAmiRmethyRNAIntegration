library(rentrez)
library(tidyverse)

#get associated genes
get_features=function(term,db){
	features=entrez_search(db=db,
					term=paste(term,"Homo sapiens[porgn]",sep=" AND "))
	features=entrez_search(db=db,
					term=paste(term,"Homo sapiens[porgn]",sep=" AND "),
					retmax=features$count)
	results=cbind(term=term,
				id=features$ids)	
return(results)}
genes=get_features("'secondary leukemia'",#search the exact phrase
				db="gene")#0 genes
genes=get_features("leukemia",db="gene")
genes=data.frame(rbind(genes,get_features("myelodysplasia",
										db="gene")))
genes%>%count(term)
#            term    n
#1       leukemia 3058
#2 myelodysplasia   39
#get gene names
library(biomaRt)
mart=useEnsembl("ensembl",
				dataset="hsapiens_gene_ensembl",
				version=106)#Ensembl Apr 2022
myannot=getBM(attributes = c("entrezgene_id","hgnc_symbol"),
			filters="entrezgene_id",
			values=genes$id,
			mart=mart)
colnames(genes)[2]="entrezgene_id"
genes=merge(genes,myannot,by="entrezgene_id",all.x=T)
genes%>%filter(!is.na(hgnc_symbol))%>%count(term)
#            term   n
#1       leukemia 2728
#2 myelodysplasia  34
genes%>%filter(!is.na(hgnc_symbol))%>%distinct(hgnc_symbol)%>%count
#     n
#1 2720
###ya son muchos, busca papers con los 2 terminos?????

########################################################
#returned ids link to https://www.ncbi.nlm.nih.gov/protein/ID
#map to entrezgene_id???????????
#protes=get_features("'secondary leukaemia'",#search the exact phrase
#				db="protein")#0 proteins
#protes=get_features("leukaemia",db="protein")
#protes=data.frame(rbind(protes,get_features("myelodysplasia",
#										db="protein")))
#protes%>%count(term)
#            term    n
#1      leukaemia 1376
#2 myelodysplasia   98
########################################################
library(clusterProfiler)

complexset=genes%>%filter(entrezgene_id!="NA")
complexset$entrezgene_id=as.character(complexset$entrezgene_id)

KEGGenrich=compareCluster(entrezgene_id~term,
	data=complexset,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
KEGGenrich=as.data.frame(KEGGenrich)

#get the list of pathways
