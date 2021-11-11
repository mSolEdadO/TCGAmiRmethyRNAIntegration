library(tidyverse)

files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)
#######################################PLOT

total=lapply(sets,function(x) x%>%count(component,omic))
total=do.call(rbind,lapply(1:5,function(x) 
	cbind(names(total)[x],total[[x]])))
colnames(total)[1]="subtype"
total$component=as.numeric(gsub("comp ","",total$component))
total$omic=factor(total$omic,levels=c("CpGs","transcripts","miRNAs"))
plots=lapply(levels(total$omic),function(x)  
	ggplot(total[total$omic==x,],aes(subtype,component,fill=n))+
	geom_tile()+xlab("")+ggtitle(x)+theme(text=element_text(size=18),
		legend.position=c(0.9,0.8),panel.background=element_blank(),
		axis.ticks.x= element_blank())+scale_fill_gradient(name="count",
		trans="log",low=1,high=400,breaks=scales::extended_breaks(n=4)))
png("../Desktop/selected.png",width=1000)
 grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=3)
dev.off()

#######################################CHECK ENRICHMENTS
sets=lapply(sets,function(x) 
	filter(x,omic=="transcripts")%>%
	select(c("variable","component")))

library(biomaRt)
library(org.Hs.eg.db)
library(KEGGREST)
library(HTSanalyzeR)
library("stringr")#do I need this???
library(GO.db)
library(GSEABase)

#get dbs 
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"entrezgene_accession"), mart=mart)
universo=as.character(myannot$entrezgene_id[!is.na(myannot$entrezgene_id)])
#[1] 29107 ids
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))#list of entrez_ids                   
#[1] 12319 BPs
name <- keggList("pathway", "hsa")
#[1] 345 pathways
names(name)=gsub("path:","",names(name))
GS_KEGG<-sapply(names(name),function(x) 
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
#61 genes have no entrezgene_id
#94 genes map 2 entrezgene_id
#$LumB
#         variable entrezgene_id entrezgene_accession
#  ENSG00000004866         93655              ST7-OT3
#  ENSG00000004866          7982                  ST7

#make diccionario... mergable
name=gsub(" - Homo sapiens (human)","",name,fixed=T)
name=as.data.frame(name)
name$hsa=rownames(name)

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
	resus=merge(resus,name,by="hsa")
})))

#enrichment of biological processes
enrichBP=lapply(sets,function(y)
	do.call(rbind,lapply(unique(y$component),function(x) {
	#lapply per comp coz I couldn't make it with group_map	
	resus=multiHyperGeoTest(
		collectionOfGeneSets=GS_GO_BP,
		universe=universo,
		hits=as.character(y$entrezgene_id[y$component==x]),
		pAdjustMethod="fdr")
	resus=cbind(resus,x)
})))
#somehow I loose multiHyper... output if these is inside the lapply
enrichBP=lapply(enrichBP,function(x) cbind(rownames(x),x))
enrichBP=lapply(enrichBP,function(x) cbind(x,Term(x[,1])))

#################CHECK OUTPUT
#pimp tables
enrichKEGG=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(enrichKEGG)[x],enrichKEGG[[x]]))))
colnames(enrichKEGG)[1]="subtype"
write_tsv(enrichKEGG,"KEGG1.enrichment")
enrichBP=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(enrichBP)[x],enrichBP[[x]]))))
colnames(enrichBP)[c(1,2,10,11)]=c("subtype","GO","component","name")
write_tsv(enrichBP,"BP1.enrichment")

enriched=list(BP=enrichBP,KEGG=enrichKEGG)
enriched=lapply(enriched,function(x) 
	x%>%filter(Adjusted.Pvalue<10**-4))

#############
#BP categories
#bps=unique(enriched$BP$GO)
#myCollection<-GOCollection(unlist(bps))
#fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
#slim <- getOBOCollection(fl)
#bpGroups=goSlim(myCollection, slim, "BP")
#% is the frequency of identifiers classified to each term

#################
#################
#plot intersections
library(igraph)
#matrix subtypes vs function
temp=lapply(enriched,function(x) x%>%distinct(name,subtype)%>%
	table%>%t)#t() so columns are
#count exclusive functions
exclusive=lapply(temp,function(y) 
 rev(table(apply(y[which(rowSums(y==0)==4),],1,paste,collapse=""))))#get all intersections
#$BP
#10000 01000 00100 00010 00001 
#   34     4    18    21    32 
#$KEGG
#10000 01000 00100 00010 00001 
#   22     1     4    12    20 
g=lapply(temp,function(x) crossprod(x))
#$BP
#        subtype
#subtype  Basal Her2 LumA LumB Normal
#  Basal     53    3    9    6     10
#  Her2       3   13    5    1      5
#  LumA       9    5   38    7     11
#  LumB       6    1    7   36      6
#  Normal    10    5   11    6     54
#
#$KEGG
#        subtype
#subtype  Basal Her2 LumA LumB Normal
#  Basal     45    4    7   13     17
#  Her2       4    9    2    2      7
#  LumA       7    2   11    4      6
#  LumB      13    2    4   26      9
#  Normal    17    7    6    9     42
values=lapply(1:2,function(y) lapply(1:5,function(x) 
	c(diag(g[[y]])[x]-exclusive[[y]][x],
		exclusive[[y]][x])))
#graph from the matrix
g1=lapply(g,function(x) 
	graph.adjacency(x,weighted=T,mode="undirected",diag=F))
V(g1$BP)$size=diag(g$BP)
V(g1$KEGG)$size=diag(g$KEGG)
pdf("enrichment.pdf")
lapply(1:2,function(x) 
	plot(g1[[x]],vertex.shape="pie",vertex.pie=values[[x]],
	vertex.size=V(g1[[x]])$size,
	edge.width=E(g1[[x]])$weight*3,
	vertex.pie.color=list(c("gray","red")),
	edge.label=E(g1[[x]])$weight,vertex.frame.color="white",
	vertex.label.cex=1.5,vertex.label.color="black",
	edge.label.cex=1.5,edge.label.color="black",main=names(g1)[x]))
dev.off()

##LO QUE SIGUE
#1) RED DE MI DE UNA FUNCION
#2)CONCLUS PROMETEDORA ¿DE LA RED?
#ENVIAR