library(tidyverse)

files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)
#######################################PLOT
library(ggplot2)
library(gridExtra)

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
		trans="log",low="#3D4E5A",high="gray",breaks=scales::extended_breaks(n=4)))
png("../Desktop/selected.png",width=1000)
 grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=3)
dev.off()

#######################################CHECK ENRICHMENTS
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

sets=do.call(rbind,lapply(1:5,function(x) cbind(names(sets)[x],sets[[x]])))
colnames(sets)[1]="subtype"
sets=sets%>%filter(omic=="transcripts")%>%
	dplyr::select(subtype,component,variable)

BPenrich=compareCluster(variable~subtype+component,
	data=sets,
	fun="enrichGO",
	OrgDb=org.Hs.eg.db,
	keyType="ENSEMBL",
	ont="BP",
	readable=T,
	pAdjustMethod = "fdr",
    pvalueCutoff  = 0.01)#slooow
write_tsv(as.data.frame(BPenrich),"BP.enrichment")

#KEGG enrichmnet needs ncbi-geneid, ncbi-proteinid or uniprot
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"entrezgene_accession"), mart=mart)
colnames(myannot)[1]="variable"
sets=merge(sets,myannot,by="variable")
	#this merge shows ensembl ids mapping to several entrez:
	#           variable entrezgene_id entrezgene_accession
	#300 ENSG00000159216           861                RUNX1
	#301 ENSG00000159216     100506403         LOC100506403
KEGGenrich=compareCluster(entrezgene_id~subtype+component,
	data=sets,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
KEGGenrich = setReadable(KEGGenrich, OrgDb = org.Hs.eg.db,
 keyType="ENTREZID")
write_tsv(as.data.frame(KEGGenrich),"KEGG.enrichment")                  

#################
#plot intersections
library(UpSetR)

enriched=list(BP=BPenrich,KEGG=KEGGenrich)
#matrix subtypes vs function
get_sets=function(enriched_table,exclusive){
 g=table(unique(enriched_table[,c("ID","subtype")]))
 if(exclusive==T){
 	g=g[rowSums(g==0)==4,]
 }
 sets=apply(g,2,function(y) names(which(y>0)))
return(sets)}
functions=lapply(enriched,get_sets,exclusive=F)
#sapply(functions,function(x) sapply(x,length))
#        BP KEGG
#Basal  471   60
#Her2   126    9
#LumA   290   12
#LumB   444   22
#Normal 616   25
pdf("enrichment.pdf")
 lapply(functions,function(x) upset(fromList(x),order.by="freq",
 	text.scale=rep(1.5,6)))
dev.off()

#check exclusive functions
#exclusive=lapply(enriched,get_sets,exclusive=T)
#shared=list()
#shared$BP=lapply(1:5,function(x) 
#	functions$BP[[x]][!functions$BP[[x]]%in%exclusive$BP[[x]]])
#shared$KEGG=lapply(1:5,function(x) 
#	functions$KEGG[[x]][!functions$KEGG[[x]]%in%exclusive$KEGG[[x]]])


##############BP categories
library(GSEABase)
library(GO.db)

# as in https://support.bioconductor.org/p/128407/
#and https://support.bioconductor.org/p/83375/
fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
#subset used for humans in PMC6800510
slim <- getOBOCollection(fl)#53 ids only
#df = select(GO.db, keys(GO.db), "ONTOLOGY")#ontology of all GOids
#table(df$ONTOLOGY[df$GOID%in%ids(slim)])#found all slim ids
#BP CC MF 
#21 16 16 
gomap=as.list(GOBPOFFSPRING)#descendents off every goid
#found=names(gomap)[names(gomap)%in%ids(slim)]
#[1] 21
#sum(found%in%df$GOID[df$ONTOLOGY=="BP"])
#[1] 21 #actually only descendents of BP goids
gomap=gomap[names(gomap)%in%ids(slim)]
slim=as.data.frame(do.call(rbind,lapply(1:21,function(x) cbind(names(gomap)[x],gomap[[x]]))))
colnames(slim)=c("parent","child")
slimnames=as.data.frame(sapply(unique(slim$PARENT),function(x) 
	Term(GOTERM[[x]])))
slimnames$id=rownames(slimnames)
colnames(slimnames)[1]="name"
slimnames=slimnames[,2:1]

BPslimenrich=compareCluster(ID~subtype+component,
	data=BPenrich,
	fun="enricher",
	TERM2GENE=slim,
	TERM2NAME=slimnames,
    pvalueCutoff  = 1)
dotplot(BPslimenrich,x="subtype",size="count")#cagada de resultado


png("BPslim.png")
ggplot(bpGroups,aes(subtype,Term,fill=Percent))+geom_tile()+
scale_fill_gradient(low="blue", high="red",na.value="white",
	trans="log10")+xlab("")+ylab("")+theme(text=element_text(size=18),
	axis.text.x = element_text(angle = 45),legend.position=c(-1,.8))
dev.off()

##LO QUE SIGUE
#1) RED DE MI DE UNA FUNCION
#2)CONCLUS PROMETEDORA ¿DE LA RED?
#ENVIAR