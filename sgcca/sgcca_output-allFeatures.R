library(tidyverse)
library(biomaRt)

#get data
files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)

#map every omic id to entrez
methy=read_tsv("../Downloads/MapMethy.tsv")
methy=methy[,c(1,4)]#keep only ids
methy=separate_rows(methy,UCSC_RefGene_Name,convert=T)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"mirbase_id","hgnc_symbol"), mart=mart)
sum(!unique(methy$UCSC_RefGene_Name)%in%myannot$hgnc_symbol)
#[1] 3549 nothing to do about them
#this aint ideal but recovers more features than entrezgene_accession
colnames(methy)[2]="hgnc_symbol"
myannot=merge(myannot,methy,by="hgnc_symbol")
#loose names, keep IDs and DB
myannot=myannot%>%dplyr::select(-1)%>%
	pivot_longer(-2,values_to="variable",names_to="ID")%>%unique()

#add entrezgene_id to selected features
sets=lapply(sets,function(x) merge(x,myannot,by="variable",all.x=T))

#########################ACTUAL ENRICHMENT
library(clusterProfiler)
library(org.Hs.eg.db)

complexset=do.call(rbind,lapply(1:5,function(x) 
	cbind(names(sets)[x],sets[[x]][,c("component","entrezgene_id","ID")])))
colnames(complexset)[1]="subtype"
#lots of features don't have an entrezgene_id
sum(is.na(complexset$entrezgene_id))
#[1] 91361
sum(!is.na(complexset$entrezgene_id))
#[1] 192023
complexset=complexset%>%filter(entrezgene_id!="NA")

BPenrich$Basal=compareCluster(entrezgene_id~component,
	data=complexset$Basal,
	fun="enrichGO",
	OrgDb=org.Hs.eg.db,
	keyType="ENSEMBL",
	ont="BP",
	readable=T,
	pAdjustMethod = "fdr",
    pvalueCutoff  = 0.01)#slooow
#No enrichment found in any of gene cluster, please check your input...
#write_tsv(as.data.frame(BPenrich),"BP-allFeatures.enrichment")

KEGGenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
KEGGenrich = as.data.frame(setReadable(KEGGenrich, OrgDb = org.Hs.eg.db,
 keyType="ENTREZID"))
write_tsv(KEGGenrich,"KEGG-allFeatures.enrichment")                  

##################PLOT INTERSECTIONS
library(UpSetR)

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
functions=get_sets(KEGGenrich,exclusive=F)
#sapply(functions,length)
# Basal   Her2   LumA   LumB Normal 
#    25     25     57     30     13 
png("enrichment-allFeatures.png")
 upset(fromList(functions),order.by="freq",
 	text.scale=rep(1.5,6))
dev.off()

#############GROUP EXCLUSIVE FUNCTIONS
exclusive=get_sets(KEGGenrich,exclusive=T)
ids=read_tsv("KEGG.exclusive")#manual classification of unique(unlist(exclusive$KEGG))
#data frame it
KEGG.classes=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(exclusive)[x],exclusive[[x]]))))
colnames(KEGG.classes)=c("subtype","ID")
KEGG.classes=merge(KEGG.classes,ids,by="ID",all.X=T)

exclusiveori=get_sets(keggori,exclusive=T)
classesori=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(exclusiveori)[x],exclusiveori[[x]]))))
colnames(classesori)=c("subtype","ID")
classesori=merge(classesori,ids,by="ID",all.X=T)
png("KEGGexclusive.png")
classesori%>%count(subtype,class)%>%
 ggplot(aes(x=n,y=class,fill=subtype))+
 geom_bar(stat="identity",position="fill")+
 annotate("text",x=1.05,y=sort(unique(classesori$class)),
 	label=classesori%>%count(class)%>%select(n)%>%unlist)+
 scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank())+xlab("")+ylab("")+
scale_fill_viridis_d(option = "plasma")
dev.off()
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

png("KEGGexclusive-allFeatures.png")
KEGG.classes%>%count(subtype,class)%>%
 ggplot(aes(x=n,y=class,fill=subtype))+
 geom_bar(stat="identity",position="fill")+
 annotate("text",x=1.05,y=sort(unique(KEGG.classes$class)),
 	label=KEGG.classes%>%count(class)%>%select(n)%>%unlist)+
 scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank())+xlab("")+ylab("")+
scale_fill_viridis_d(option = "plasma")
dev.off()

#############ADD GSEA INFO
library(ggrepel)

gsea=read_tsv("KEGG.gsea")
gsea=gsea%>%select(subtype,Description,NES,p.adjust)
gsea=cbind(gsea[,2:4],subtype=gsub("Luma","LumA",gsub("Lumb","LumB",
		str_to_title(gsub("_normal","",gsea$subtype)))))
#count enriched components per function
temp=KEGGenrich%>%group_by(subtype,Description)%>%tally
temp=merge(temp,gsea,by=c("subtype","Description"))
png("NES-ncomp-allFeatures.png")
ggplot(temp,aes(y=NES,x=n,alpha=-log(p.adjust),color=subtype))+
	geom_point(size=3)+xlab("components")+theme_light(base_size=18)+
	scale_color_manual(values=c("#0D0887","#7E03A8","#CC4678","#F89441"))+
	scale_alpha(range = c(0.1,1),breaks=c(0,3,15))+
	scale_x_continuous(breaks=seq(0,max(temp$n),2))+
	geom_text_repel(aes(label=ifelse((p.adjust<0.001)&(n>4),
	Description,'')),alpha=1)
dev.off()
j=temp$Description[temp$n>1&temp$p.adjust<0.05]
#30 

#############CHECK SHARED FUNCTIONS
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
png("KEGGenrichment-allFeatures.png",width=800,height=600)
ggplot(KEGGshared)+geom_point(aes(x=subtype,y=Description,
	size=components,col=genes))+xlab("")+ylab("")+
	scale_color_gradient(low="blue",high="red")+
	theme_light(base_size=18)+scale_size(range=c(2,10))+
	theme(axis.ticks=element_blank())+
#	annotate("text",y=j[j%in%KEGGshared$Description],x="Basal",
#		label="*",size=7,vjust=.8,hjust=5)+
	annotate("text",y=j[j%in%KEGGshared$Description],x="Normal",
		label="*",size=7,vjust=.8,hjust=-4.1)+
	coord_cartesian(clip="off")
dev.off()

