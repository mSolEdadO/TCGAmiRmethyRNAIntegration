library(tidyverse)
library(biomaRt)

#get data
files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)

#map every omic id to entrez
methy=read_tsv("../Downloads/MapMethy.tsv")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"mirbase_id"), mart=mart)
#prepare to merge with sets
myannot=myannot%>%pivot_longer(-entrezgene_id,names_to="DB",
	values_to="variable")
myannot=myannot%>%dplyr::select(entrezgene_id,variable)%>%
	filter(variable!="")
methy=methy%>%dplyr::select(entrezgene_id,IlmnID)
colnames(methy)[2]="variable"
myannot=rbind(myannot,methy)
#add entrezgene_id to selected features
sets=lapply(sets,function(x) merge(x,myannot,all.x=T,by="variable"))
#all the transcripts serve for enrichment, some miRNAs & CpGs don't
#sapply(setsAlt,function(x) table(x$omic[is.na(x$entrezgene_id)]))
#       Basal Her2  LumA  LumB Normal
#CpGs   14123 4038 41035 13250  10200
#miRNAs     6    1    25     7      3

#########################ACTUAL ENRICHMENT
library(clusterProfiler)
library(org.Hs.eg.db)

complexset=do.call(rbind,lapply(1:5,function(x) 
	cbind(names(sets)[x],sets[[x]][,c("component","entrezgene_id","omic")])))
colnames(complexset)[1]="subtype"
complexset=complexset%>%filter(entrezgene_id!="NA")
complexset$entrezgene_id=as.character(complexset$entrezgene_id)

BPenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichGO",
	OrgDb=org.Hs.eg.db,
	ont="BP",
	readable=T,
	pAdjustMethod = "fdr",
    pvalueCutoff  = 0.01)#slooow
write_tsv(as.data.frame(BPenrich),"BP-allFeatures.enrichment")

KEGGenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
write_tsv(as.data.frame(KEGGenrich),"KEGG-allFeatures.enrichment")                  

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
#Basal  226   26
#Her2   127   31
#LumA   820   72
#LumB   448   35
#Normal  48   26
pdf("enrichment-allFeatures.pdf")
 lapply(functions,function(x) upset(fromList(x),order.by="freq",
 	text.scale=rep(1.5,6)))
dev.off()

#############GROUP EXCLUSIVE FUNCTIONS
exclusive=lapply(enriched,get_sets,exclusive=T)
#write_tsv(unique(KEGGenrich[KEGGenrich$ID%in%unlist(exclusive$KEGG),c("ID","Description")]),"KEGG.exclusive")
#class comes from https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=
#kegg=readLines("br08901.keg")
#kegg=kegg[grep("^[A-Z]",kegg,perl=T)]
#kegg=kegg[substr(kegg,1,1)!="A"]
#i=which(substr(kegg,1,1)=="B")
#i=c(i,length(kegg)+1)
#keggL=lapply(1:59,function(x) kegg[(i[x]+1):(i[x+1]-1)])
#names(keggL)=substr(kegg[i[1:59]],4,nchar(kegg[i[1:59]]))
#keggL=lapply(keggL,function(x) substr(x,13,nchar(x)))
#keggDF=as.data.frame(do.call(rbind,lapply(1:length(keggL),function(x) cbind(names(keggL)[x],keggL[[x]]))))
#colnames(keggDF)=c("class","Description")
#merge("KEGG.exclusive",keggDF,by="Description",all.x=T)
#also simplified classes with ':' to the first term
ids=read_tsv("KEGG.exclusive")
#data frame it
KEGG.classes=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(exclusive$KEGG)[x],exclusive$KEGG[[x]]))))
colnames(KEGG.classes)=c("subtype","ID")
KEGG.classes=merge(KEGG.classes,ids,by="ID")

#get BP categories
library(GSEABase)
library(GO.db)
# as in https://support.bioconductor.org/p/128407/
#and https://support.bioconductor.org/p/83375/
fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
#subset used for humans in PMC6800510
slim <- getOBOCollection(fl)#53 ids only
df = select(GO.db, keys(GO.db), "ONTOLOGY")#ontology of all GOids
table(df$ONTOLOGY[df$GOID%in%ids(slim)])#found all slim ids
#BP CC MF 
#21 16 16 
gomap=as.list(GOBPOFFSPRING)#descendents off every goid
found=names(gomap)[names(gomap)%in%ids(slim)]
#[1] 21
sum(found%in%df$GOID[df$ONTOLOGY=="BP"])
#[1] 21 #actually only descendents of BP goids
gomap=gomap[names(gomap)%in%ids(slim)]
#format to easy data frames
slim=as.data.frame(do.call(rbind,lapply(1:21,function(x) 
	cbind(names(gomap)[x],gomap[[x]]))))
colnames(slim)=c("parent","child")
slimnames=as.data.frame(sapply(unique(slim$parent),function(x) 
	Term(GOTERM[[x]])))
slimnames$parent=rownames(slimnames)
colnames(slimnames)[1]="name"
#count categories per subtype
BP.classes=as.data.frame(do.call(rbind,lapply(1:5,
	function(x) cbind(names(exclusive$BP)[x],exclusive$BP[[x]]))))
colnames(BP.classes)=c("subtype","child")
BP.classes=merge(merge(BP.classes,slim,by="child"),
	slimnames,by="parent")

############TEST OVER-REPRESENTATION OF EXCLUSIVE FUNCTIONS
library(ggplot2)
bias=function(classes,subtype,class){
	totals=table(classes[,c(class,subtype)])
	s=colSums(totals)
	ps=p.adjust(apply(totals,1,function(x) 
	fisher.test(rbind(x,s-x),simulate.p.value=T)$p.val))
	#TRUE needed coz > 2 categories
	return(rownames(totals)[ps<0.05])}
bias(KEGG.classes,2,4)
#character(0)
i=bias(BP.classes,3,4)
#[1] "DNA metabolic process"   "immune system process"  
#[3] "lipid metabolic process"

png("KEGGexclusive-allFeatures.png",width=600)
KEGG.classes%>%count(subtype,class)%>%
 ggplot(aes(x=n,y=class,fill=subtype))+
 geom_bar(stat="identity",position="fill")+
 annotate("text",x=1.05,y=sort(unique(KEGG.classes$class)),
 	label=KEGG.classes%>%count(class)%>%dplyr::select(n)%>%unlist)+
 scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank())+xlab("")+ylab("")+
scale_fill_viridis_d(option = "plasma")
dev.off()
png("BPexclusive-allFeatures.png",width=600)
BP.classes%>%count(subtype,name)%>%
ggplot(aes(x=n,y=name,fill=subtype))+
geom_bar(stat="identity",position="fill")+
scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank())+xlab("")+ylab("")+
scale_fill_viridis_d(option = "plasma")+
annotate("text",x=1.05,y=sort(unique(BP.classes$name)),
 	label=BP.classes%>%count(name)%>%dplyr::select(n)%>%unlist)+
annotate("text",y=i,x=-.05,label="*",size=8,vjust=.8)
dev.off()

#############ADD GSEA INFO
#library(ggrepel)

files=list.files()
files=files[grep("gsea",files)]
gsea=lapply(files,read_tsv)
names(gsea)=gsub(".gsea","",files)
#match gsea subtype to SGCCA result
gsea=lapply(gsea,function(x) 
	x%>%dplyr::select(subtype,Description,NES,p.adjust))
gsea=lapply(gsea,function(x) 
	cbind(subtype=gsub("Luma","LumA",gsub("Lumb","LumB",
		str_to_title(gsub("_normal","",x$subtype)))),x[,2:4]))
#count enriched components per function
temp=lapply(enriched,function(x)
	x%>%group_by(subtype,Description)%>%tally)
temp=lapply(1:2,function(x) merge(temp[[x]],gsea[[x]],
	by=c("subtype","Description")))

pdf("NES-ncomp-allFeatures.pdf")
lapply(temp,function(x)
	ggplot(x,aes(y=NES,x=n,alpha=-log(p.adjust),color=subtype))+
	geom_point(size=3)+xlab("components")+theme_light(base_size=18)+
	scale_color_manual(values=c("#0D0887","#7E03A8","#CC4678","#F89441"))+
	scale_alpha(range = c(0.1,1),breaks=c(0,3,15))+
	geom_text_repel(aes(label=ifelse((p.adjust<0.001)&(n>(max(n)/3)),
	Description,'')),alpha=1))
dev.off()
j=lapply(temp,function(x)
	unique(x$Description[x$n>1&x$p.adjust<0.05]))
#  BP KEGG 
# 233   28

#############CHECK SHARED FUNCTIONS
heatmatrix=function(enrichment){
	#only functions enriched in more than one dataset
	i=enrichment%>%distinct(subtype,Description)%>%count(Description)%>%
 	filter(n>1)%>%dplyr::select(Description)%>%unlist
	enrichment=enrichment[enrichment$Description%in%i,]
	#count genes & component per function and subtype
	edges=enrichment%>%dplyr::select(subtype,Description,geneID)%>%
		group_by(subtype,Description)
	edges=edges%>%group_map(~length(unique(unlist(strsplit(.x$geneID,
		"/")))))%>%unlist%>%cbind(group_keys(edges))
	colnames(edges)[1]="genes"
	edges1=enrichment%>%group_by(subtype,Description)%>%tally
	colnames(edges1)[3]="components"
	edges=merge(edges,edges1,by=c("subtype","Description"))
return(edges)}
shared=lapply(enriched,heatmatrix)
#sapply(shared,function(x) length(unique(x$Description)))
#  BP KEGG 
# 403   58 
#BP is too large to plot

png("KEGGenrichment-allFeatures.png",width=800,height=700)
ggplot(shared$KEGG)+geom_point(aes(x=subtype,y=Description,
	size=components,col=genes))+xlab("")+ylab("")+
	scale_color_gradient(low="blue",high="red")+
	theme_light(base_size=18)+scale_size(range=c(2,10))+
	theme(axis.ticks=element_blank())+
	annotate("text",y=j$KEGG[j$KEGG%in%shared$KEGG$Description],x="Basal",
		label="*",size=7,vjust=.8,hjust=4.5)+
	coord_cartesian(clip="off")
dev.off()
#plot the 19 BPs found in the 5 datasets
temp=names(which(table(shared$BP$Description)==5))
#[1] 19
png("BPenrichment-allFeatures.png",width=850)
ggplot(shared$BP[shared$BP$Description%in%temp,])+geom_point(aes(x=subtype,y=Description,
	size=components,col=genes))+xlab("")+ylab("")+
	scale_color_gradient(low="blue",high="red")+
	theme_light(base_size=18)+scale_size(range=c(2,10))+
	theme(axis.ticks=element_blank())
dev.off()
#sum(j$BP%in%temp)
#[1] 19 all affected by differential expression
#############COMPARE WITH TRANSCRIPTS ONLY ENRICHMENT
files=list.files()
files=files[grep(".enrichment",files,fixed=T)]
files=files[grep("all",files,invert=T)]
ori=lapply(files,read_tsv)
names(ori)=gsub(".enrichment","",files)
temp=as.data.frame(cbind(set="transcripts only",
	p.adjust=unlist(sapply(ori,function(x) 
		x%>%dplyr::select(p.adjust)))))
temp=rbind(temp,as.data.frame(cbind(set="all Features",
	p.adjust=unlist(sapply(enriched,function(x) 
		x%>%dplyr::select(p.adjust))))))
temp$p.adjust=as.numeric(temp$p.adjust)
temp$DB=sapply(strsplit(rownames(temp),".",fixed=T),function(x) x[1])
png("p.values.png")
ggplot(temp,aes(x=p.adjust,col=set,fill=set))+
geom_density(aes(y=..scaled..),alpha=0.3)+
facet_wrap(~DB)+theme(text=element_text(size=18))+
scale_x_continuous(breaks=c(0,0.005,0.01))
dev.off()