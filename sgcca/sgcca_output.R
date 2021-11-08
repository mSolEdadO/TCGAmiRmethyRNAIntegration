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
plots=lapply(unique(total$omic),function(x)  
	ggplot(total[total$omic==x,],aes(subtype,component,fill=n))+
	geom_tile()+xlab("")+ggtitle(x)+theme(text=element_text(size=18))+
	scale_fill_gradient(name = "count",trans="log",low=1,high=400))
#still wanna:
#1) move color key to the bottom,
#2) erase background grid
#3) 5 breaks in color key, no matter omic
grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=3)



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
enrichKEGG=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(enrichKEGG)[x],enrichKEGG[[x]]))))
colnames(enrichKEGG)[1]="subtype"
write_tsv(enrichKEGG,"KEGG1.enrichment")
temp=enrichKEGG%>%filter(Adjusted.Pvalue<0.01)%>%
	distinct(subtype,hsa)
#pathways per subtype
temp%>%count(subtype)
# subtype   n
#   Basal 269
#    Her2 204 
#    LumA 229
#    LumB 264
#  Normal 249
paths=temp%>%count(hsa)%>%filter(n==5)%>%dplyr::select(hsa)
# 165 pathways enriched for all the datasets
#enrichKEGG%>%filter(hsa%in%temp$hsa)%>%dplyr::select(name)%>%distinct()
targets=myannot$variable[myannot$entrezgene_id%in%
	unlist(GS_KEGG_id[names(GS_KEGG_id)%in%paths$hsa])]
#Ensembl id of -2541- genes linked to the 55 common pathways

enrichBP=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(enrichBP)[x],enrichBP[[x]]))))
colnames(enrichBP)[c(1,2,10,11)]=c("subtype","GO","component","name")
write_tsv(enrichBP,"BP1.enrichment")
temp=enrichBP%>%filter(Adjusted.Pvalue<0.01)%>%
	distinct(subtype,GO)
temp%>%count(subtype)
# subtype    n
#   Basal 1213
#    Her2  806
#    LumA  958
#    LumB 1255
#  Normal 1047
bps=temp%>%count(GO)%>%filter(n==5)%>%dplyr::select(GO)
# 496
targets=union(targets,myannot$variable[myannot$entrezgene_id%in%
	unlist(GS_GO_BP[names(GS_GO_BP)%in%bps$GO])])
#add the 4867 ensembl ids linked to the 138 common bps
#this are the targets for MI 

#BP categories
myCollection<-GOCollection(unlist(bps))
fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
slim <- getOBOCollection(fl)
bpGroups=goSlim(myCollection, slim, "BP")
#% is the frequency of identifiers classified to each term


#########################plots to refine
#png("selected.png")
#heatmap.2(A1[rev(rownames(A)),],scale='n',trace='n',dendrogram='n',Colv=F,Rowv=F,col=rev(grey.colors(20)),key=T,labRow="",srtCol=0,labCol=c("","Basal","","","Her2","","","LumA","","","LumB","","","Normal",""),lmat=rbind(c(0,3),c(2,1),c(0,4)),lhei=c(0.5,5,1.3),lwid=c(0.2,8),colsep=c(3,6,9,12),key.title="",cexCol=2,density.info='n',key.xlab="features",margins=c(5,6),adjCol=c(0.5,NA),add.expr=text(labels=apply(A[,2:ncol(A)],2,function(x) paste(min(as.numeric(x),na.rm=T),max(as.numeric(x),na.rm=T),sep=' - ')),y= as.numeric(sapply(c(36,21,60,40,19),rep,3)),srt=45,xpd=NA,x=(1:15)+.18,cex=1.2))
#axis(side=4,at=seq(.19,1,.15),labels=seq(0,50,10),line=-1)
#mtext(side=4,"components",line=1)
#dev.off()
