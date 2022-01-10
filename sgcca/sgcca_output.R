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
total$component=as.numeric(gsub("comp*","",total$component,perl=T))
total$omic=factor(total$omic,levels=c("CpGs","transcripts","miRNAs"))
plots=lapply(levels(total$omic),function(x)  
	ggplot(total[total$omic==x,],aes(subtype,component,fill=n))+
	geom_tile()+xlab("")+ggtitle(x)+theme(text=element_text(size=18),
		legend.position=c(0.9,0.8),panel.background=element_blank(),
		axis.ticks.x= element_blank(),
		axis.text.x=element_text(vjust=5))+
		scale_fill_viridis_b(name="count"))
png("selected.png",width=880,height=400)
 grid.arrange(plots[[1]],plots[[2]],plots[[3]],ncol=3)
dev.off()

#######################################ENRICHMENT
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(enrichplot)

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

##################PLOT INTERSECTIONS
library(UpSetR)

enriched=list(BP=BPenrich,KEGG=KEGGenrich)
get_sets=function(enriched_table,exclusive){
#matrix subtypes vs function
 g=table(unique(enriched_table[,c("ID","subtype")]))
#to get exclusive functions for another time
 if(exclusive==T){
 	g=g[rowSums(g==0)==4,]
 }
 #upset() needs a list of IDs
 sets=apply(g,2,function(y) names(which(y>0)))
return(sets)}
functions=lapply(enriched,get_sets,exclusive=F)
#sapply(functions,function(x) sapply(x,length))
#        BP KEGG
#Basal  277   40
#Her2    73    8
#LumA   551   76
#LumB   247   14
#Normal 313   22
pdf("enrichment.pdf")
 lapply(functions,function(x) upset(fromList(x),order.by="freq",
 	text.scale=rep(1.5,6)))
dev.off()

#############CHECK SHARED FUNCTIONS
library(enrichplot)
heatmatrix=function(enrichment){
	#only functions enriched in more than one dataset
	i=enrichment%>%distinct(subtype,Description)%>%count(Description)%>%
 	filter(n>1)%>%select(Description)%>%unlist
	enrichment=enrichment[enrichment$Description%in%i,]
	#count genes & component per function and subtype
	edges=enrichment%>%select(subtype,Description,geneID)%>%
		group_by(subtype,Description)
	edges=edges%>%group_map(~length(unique(unlist(strsplit(.x$geneID,
		"/")))))%>%unlist%>%cbind(group_keys(edges))
	colnames(edges)[1]="genes"
	edges1=enrichment%>%group_by(subtype,Description)%>%tally
	colnames(edges1)[3]="components"
	edges=merge(edges,edges1,by=c("subtype","Description"))
return(edges)}
KEGGshared=heatmatrix(KEGGenrich)
#[1] 35 functions
png("KEGGenrichment.png")
ggplot(KEGGshared)+geom_point(aes(x=subtype,y=Description,
	size=components,col=genes))+xlab("")+
	theme(text=element_text(size=18))+
	scale_color_gradient(low="blue",high="red")+theme_light()+
	scale_size(range=c(2,10))
dev.off()

#select functions recursively enriched for further analyses
temp=KEGGshared%>%select(-genes)%>%pivot_wider(names_from="subtype",
	values_from="components")
KEGGshared=temp$Description[rowSums(as.data.frame(temp)[,2:6]>1,na.rm=T)>0]
#[1] 20 functions selected
BPshared=heatmatrix(BPenrich)#is too large, so I don't plot it
#[1] 271 functions
temp=BPshared%>%select(-genes)%>%pivot_wider(names_from="subtype",
	values_from="components")
BPshared=temp$Description[rowSums(as.data.frame(temp)[,2:6]>1,na.rm=T)>0]
#[1] 123 functions selected


#############CHECK EXCLUSIVE FUNCTIONS
#exclusive=lapply(enriched,get_sets,exclusive=T)
#manual classification of unlist(unique(exclusive$KEGG))
KEGG.classes=as.data.frame(cbind(unique(ids$class),
	do.call(rbind,lapply(unique(ids$class),function(x)
	 sapply(exclusive$KEGG,function(y) sum(y%in%ids$ID[ids$class==x]))))))
KEGG.classes=KEGG.classes%>%pivot_longer(-1,
	names_to="subtype",values_to="pathways")
KEGG.classes$pathways=as.numeric(KEGG.classes$pathways)
KEGG.classes=KEGG.classes%>%group_by(V1)%>%
	summarise(total=sum(pathways))%>%merge(KEGG.classes,by="V1")
png("KEGGexclusive.png")
 ggplot(KEGG.classes,aes(x=pathways,y=V1,
 	fill=subtype))+geom_bar(stat="identity",position="fill")+
 annotate("text",x=1.05,y=unique(KEGG.classes$V1),
 	label=KEGG.classes$total[seq(1,nrow(KEGG.classes),5)])+
 scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank())+xlab("")+ylab("")+
scale_fill_viridis_d(option = "plasma")
dev.off()

#BP categories
#library(GSEABase)
#library(GO.db)
# as in https://support.bioconductor.org/p/128407/
#and https://support.bioconductor.org/p/83375/
#fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
#subset used for humans in PMC6800510
#slim <- getOBOCollection(fl)#53 ids only
#df = select(GO.db, keys(GO.db), "ONTOLOGY")#ontology of all GOids
#table(df$ONTOLOGY[df$GOID%in%ids(slim)])#found all slim ids
#BP CC MF 
#21 16 16 
#gomap=as.list(GOBPOFFSPRING)#descendents off every goid
#found=names(gomap)[names(gomap)%in%ids(slim)]
#[1] 21
#sum(found%in%df$GOID[df$ONTOLOGY=="BP"])
#[1] 21 #actually only descendents of BP goids
#gomap=gomap[names(gomap)%in%ids(slim)]
#slim=as.data.frame(do.call(rbind,lapply(1:21,function(x) cbind(names(gomap)[x],gomap[[x]]))))
#colnames(slim)=c("parent","child")
#slimnames=as.data.frame(sapply(unique(slim$PARENT),function(x) 
#	Term(GOTERM[[x]])))
#slimnames$id=rownames(slimnames)
#colnames(slimnames)[1]="name"
#slimnames=slimnames[,2:1]
#BPslimenrich=compareCluster(ID~subtype+component,
#	data=as.data.frame(BPenrich),
#	fun="enricher",
#	TERM2GENE=slim,
#	TERM2NAME=slimnames,
#   pvalueCutoff  = 1)
#dotplot(BPslimenrich,x="subtype",size="count")#cagada de resultado
#png("BPslim.png")
#ggplot(bpGroups,aes(subtype,Term,fill=Percent))+geom_tile()+
#scale_fill_gradient(low="blue", high="red",na.value="white",
#	trans="log10")+xlab("")+ylab("")+theme(text=element_text(size=18),
#	axis.text.x = element_text(angle = 45),legend.position=c(-1,.8))
#dev.off()


#BPsem=pairwise_termsim(BPenrich,semData=godata('org.Hs.eg.db',
# ont="BP"))
#BPsem@termsim has the semantic similarity bewteen go terms
#group=treeplot(BPsem,nCluster=50,showCategory=temp)
#nCluster chosen after lots of plots
#groups=as.data.frame(cbind(group$data$label,group$data$group))
#colnames(groups)=c("Description","group")
#groups=groups[!is.na(groups$Description),]
#heatBP=merge(heatBP,groups,by="Description")
#heatBP1=heatBP%>%group_by(subtype,group)%>%summarise(genes=sum(genes),
#	components=sum(components),processess=length(unique(Description)))


#temp=clusterProfiler::simplify(BPenrich,
#	cutoff=0.01,#semantic similarity higher than `cutoff` are redundant 
#	by="p.adjust",
#	select_fun=min)#select representative term by min p.adjust?




#→→→→→→→→→→→→reduce por NES/p.value
#cuadrantes size vs NES, pvalue = transparencia, color por subtipo?
 