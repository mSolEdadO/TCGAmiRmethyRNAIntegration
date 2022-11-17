library(tidyverse)
library(UpSetR)

bp=read_tsv("BP.enrichment",show_col_types = FALSE)
k=read_tsv("KEGG.enrichment",show_col_types = FALSE)

enriched=list(BP=bp,
			  KEGG=k)
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
	fisher.test(rbind(x,s-x),simulate.p.value=T)$p.val),"bonferroni")
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
 	legend.position="bottom",legend.margin=margin(-25,0,0,-150))+
 xlab("")+ylab("")+scale_fill_viridis_d(option = "plasma")+
annotate("text",x=1.08,y=sort(unique(BP.classes$name)),
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
