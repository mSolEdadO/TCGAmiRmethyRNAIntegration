library(tidyverse)
library(biomaRt)

#get data
files=list.files()
files=files[grep("stable",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".stable","",files)

#map every omic id to entrez
methy=read_tsv("../MapMethy.tsv")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=105)
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

#########################ACTUAL ENRICHMENT#########################
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
write_tsv(as.data.frame(BPenrich),"BP.enrichment")

KEGGenrich=compareCluster(entrezgene_id~subtype+component,
	data=complexset,
	fun="enrichKEGG",
	pAdjustMethod = "fdr",
	organism = 'hsa',
	pvalueCutoff = 0.01)
#get a nice table
write_tsv(as.data.frame(KEGGenrich),"KEGG.enrichment")                  

##################PLOT INTERSECTIONS
library(UpSetR)

enriched=list(BP=as.data.frame(BPenrich),
			  KEGG=as.data.frame(KEGGenrich))
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
#Basal  105   12
#Her2    91    4
#LumA   553   54
#LumB   248   15
#Normal  44    2
pdf("enrichment.pdf",height=5)
 lapply(functions,function(x) upset(fromList(x),text.scale=1.5,order.by="degree"))
dev.off()

###################EXCLUSIVE FUNCTIONS###########3#################
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
KEGG.classes=as.data.frame(do.call(rbind,lapply(1:4,function(x) 
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
#[1] "cellular component organization" "DNA metabolic process"          
#[3] "establishment of localization"  

#ni vale la pena dibujarlo
#png("KEGGexclusive.png")
#KEGG.classes%>%count(subtype,class)%>%
# ggplot(aes(x=n,y=class,fill=subtype))+
# geom_bar(stat="identity",position="fill")+
# annotate("text",x=1.05,y=sort(unique(KEGG.classes$class)),
# 	label=KEGG.classes%>%count(class)%>%dplyr::select(n)%>%unlist)+
# scale_x_continuous(labels=scales::percent)+
# theme(text=element_text(size=18),axis.ticks=element_blank(),
# 	panel.background=element_blank(),legend.title=element_blank())+
# xlab("")+ylab("")+scale_fill_viridis_d(option = "plasma")
#dev.off()

png("BPexclusive.png")
BP.classes%>%count(subtype,name)%>%
ggplot(aes(x=n,y=name,fill=subtype))+
geom_bar(stat="identity",position="fill")+
scale_x_continuous(labels=scales::percent)+
 theme(text=element_text(size=18),axis.ticks=element_blank(),
 	panel.background=element_blank(),legend.title=element_blank(),
 	legend.position="bottom",legend.margin=margin(-20,0,0,0))+
 xlab("")+ylab("")+scale_fill_viridis_d(option = "plasma")+
annotate("text",x=1.05,y=sort(unique(BP.classes$name)),
 	label=BP.classes%>%count(name)%>%dplyr::select(n)%>%unlist)+
annotate("text",y=i,x=-.05,label="*",size=8,vjust=.8)
dev.off()

#to check exclusive categories
#temp=sapply(i,function(y) 
#	p.adjust(sapply(unique(BP.classes$subtype),function(x) fisher.test(
#		matrix(c(sum(BP.classes$name==y&BP.classes$subtype==x),
#				sum(BP.classes$name==y&BP.classes$subtype!=x),
#				sum(BP.classes$name!=y&BP.classes$subtype==x),
#				sum(BP.classes$name!=y&BP.classes$subtype!=x)),
#		ncol=2),alternative="g")$p.val)))
#       cellular component organization DNA metabolic process
#LumA                        1.00000000          1.0000000000
#Her2                        1.00000000          1.0000000000
#Basal                       0.04023048          1.0000000000
#LumB                        1.00000000          1.0000000000
#Normal                      0.04124813          0.0005866016
#       establishment of localization
#LumA                    8.607994e-06
#Her2                    1.000000e+00
#Basal                   1.000000e+00
#LumB                    1.000000e+00
#Normal                  1.000000e+00

#####################CHECK SHARED FUNCTIONS###################3

#ADD GSEA INFO
files=list.files()
files=files[grep("gsea",files)]
gsea=lapply(files,read_tsv)
names(gsea)=gsub(".gsea","",files)
#match gsea subtype to SGCCA result
gsea=lapply(gsea,function(x) 
	x%>%dplyr::select(subtype,Description,NES,p.adjust))
gsea$KEGG$subtype=gsub("_Normal","",gsea$KEGG$subtype)

j=lapply(gsea,function(x)
	unique(x$Description[x$p.adjust<0.05]))
#names(j)=names(enriched) check j is named
#  BP KEGG 
# 130   11 

#get shared function
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
# 233   16 
#BP is too large to plot

############JACCARD MATRIX
#matrix with the number of components where each pair of
# functions are enriched together	#
subtypes=unique(enriched$KEGG$subtype)
coenriched=lapply(enriched,function(x) lapply(subtypes,
	function(y) x[x$subtype==y,]))
coenriched=lapply(coenriched,function(z) lapply(z,function(w) 
	1-sapply(unique(w$Description),function(x) 
		sapply(unique(w$Description),function(y) 
			length(intersect(w$component[w$Description==x],
				w$component[w$Description==y]))/
			length(union(w$component[w$Description==x],
				w$component[w$Description==y]))))))
#u coul also just show the trees
trees=lapply(coenriched,function(x) lapply(x,function(y) 
	hclust(as.dist(y))))

#get the groups that are enriched exactly in the same components
groups=lapply(trees,function(x) lapply(x,function(y) 			
	cutree(y,h=0)))	
#paste with heatmatrix output
groups=lapply(groups,function(x) lapply(1:5,function(y) 
	data.frame(cbind("subtype"=subtypes[y],
					 "Description"=names(x[[y]]),
					 "group"=x[[y]]))))
groups=lapply(groups,function(x) do.call(rbind,x))

#for the plot
temp=lapply(names(shared),function(x) 
			merge(groups[[x]],shared[[x]],
				by=c("subtype","Description"),
				all.y=T))
#ungrouped functions should be the same figure
temp=lapply(temp,function(x) x%>%group_by(subtype)%>%
									add_count(group))
names(temp)=names(shared)
temp$BP$group[temp$BP$n==1]=0
temp$KEGG$group[temp$KEGG$n==1]=0

png("KEGGshared.png",width=600)
 ggplot(temp$KEGG)+geom_point(aes(x=subtype,y=Description,
 	size=components,col=genes,shape=group))+
 scale_color_gradient(low="blue",high="red")+
 theme_light(base_size=16)+scale_size(range=c(3,10))+
 theme(axis.ticks=element_blank())+annotate("text",
 	y=j$KEGG[j$KEGG%in%shared$KEGG$Description],x="Basal",label="*",
 	size=7,vjust=.8,hjust=4.3)+coord_cartesian(clip="off")+
 scale_shape_manual(values=c(16:18,15),guide="none")+xlab("")+ylab("")
dev.off()

#just draw bps found >= 3 datasets or it'll be too large 
temp$BP=temp$BP%>%filter(Description%in%names(which(rowSums(table(temp[[1]][,2:1]))>3)))
#trick to repeat figures between datasets
temp$BP$group=temp$BP%>%group_map(~(0:8)[as.factor(.x$group)])%>%unlist
png("BPeshared.png",width=850,height=700)
 ggplot(temp$BP)+geom_point(aes(x=subtype,y=Description,
 	size=components,col=genes,shape=as.character(group)))+
 scale_color_gradient(low="blue",high="red")+
 theme_light(base_size=16)+scale_size(range=c(3,10),trans="log10")+
 theme(axis.ticks=element_blank())+xlab("")+ylab("")+
 scale_shape_manual(values=c(16:18,15,0:2,5,6),guide="none")+
 annotate("text",y=j$BP[j$BP%in%temp$BP$Description],x="Basal",
 	label="*",size=7,vjust=.8,hjust=5.3)+coord_cartesian(clip="off")
dev.off()
