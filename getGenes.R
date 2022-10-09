library(rentrez)
library(tidyverse)

#searching directly on gene db has unspecific results
get_features=function(term){
	#get related papers
	papers=entrez_search(db="pubmed",
	term=term)
	papers=entrez_search(db="pubmed",
	term=term,
	retmax=papers$count)
	#get the genes from those papers
	ids <- entrez_link(dbfrom='pubmed', 
		id=papers$ids,
		db='gene')
	results=cbind(term=term,
				entrezgene_id=ids$links$pubmed_gene)	
return(results)}
genes=get_features('"secondary leukemia" AND MDS')
genes=data.frame(rbind(genes,
					get_features('"secondary leukemia" AND myelodysplasia')))

genes=genes%>%filter(!duplicated(entrezgene_id))
genes%>%count(term)
#                          term n
#1 "secondary leukemia" AND MDS 8
#"secondary leukemia" AND myelodysplasia outputs redundant genes

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
library(org.Hs.eg.db)

#complexset=genes%>%filter(entrezgene_id!="NA")
genes$entrezgene_id=as.character(complexset$entrezgene_id)

#simultaneous enrichment of several sets
KEGGenrich=compareCluster(entrezgene_id~term,
	data=genes,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#No enrichment found in any of gene cluster, please check your input...
BPenrich=compareCluster(entrezgene_id~term,
	data=genes,
	fun="enrichGO",
	OrgDb=org.Hs.eg.db,
	ont="BP",
	pAdjustMethod = "fdr",
    pvalueCutoff  = 0.01)#slooow
#get a nice table
BPenrich=as.data.frame(BPenrich)
# dim(BPenrich)
#[1] 29 11

map=as.data.frame(org.Hs.egGO)
colnames(BPenrich)=gsub("ID","go_id",colnames(BPenrich))
#somehow not all bps are recoverable
BPenrich$Description[!BPenrich$ID%in%map$go_id]
#[1] "response to acid chemical"                    
#[2] "formation of extrachromosomal circular DNA"   
#[3] "negative regulation of organelle organization"
#[4] "cellular response to chemical stress"         
added_genes=BPenrich%>%merge(map,by="go_id",all.x=T)%>%
			dplyr::select(gene_id)

genes=rbind(genes,cbind(entrezgene_id=unlist(added_genes),
						term="enrichment"))

#get gene symbols
library(biomaRt)
mart=useEnsembl("ensembl",
				dataset="hsapiens_gene_ensembl",
				version=106)#Ensembl Apr 2022
myannot=getBM(attributes = c("entrezgene_id","hgnc_symbol"),
			filters="entrezgene_id",
			values=genes$entrezgene_id,
			mart=mart)
genes=genes%>%filter(!is.na(entrezgene_id))%>%
			merge(myannot,by="entrezgene_id")%>%
			distinct
genes%>%filter(!is.na(hgnc_symbol))%>%count(term)
#                          term   n
#1 "secondary leukemia" AND MDS   8
#2                   enrichment 734
