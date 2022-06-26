library(tidyverse)
library(biomaRt)
library(rentrez)

####################PER SUBTYPE ANALYSIS################### 
files=list.files()
files=files[grep("selected",files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)

#up-set plot per omic & subtype of selected features
omics=unique(features$Basal$omic)
features=lapply(omics,function(x) lapply(sets,function(y)
			#list features selected
			y%>%filter(omic==x)%>%distinct(variable)%>%unlist))
pdf("slctd_feats_intrsctn.pdf")
 lapply(features,function(x) upset(fromList(x),order.by="freq"))
dev.off()

#top selected features
top=lapply(sets,function(x) x%>%count(variable,omic)
	%>%arrange(desc(n))%>%group_by(omic)%>%group_map(~head(.x,5)))
top=lapply(top,function(x) do.call(rbind,x))
top=as.data.frame(do.call(rbind,lapply(1:5,function(x)
						 cbind(top[[x]],subtype=names(top)[x]))))
#add gene descriptions
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	version=105)
#https://dec2021.archive.ensembl.org
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol",
	"wikigene_description"),mart=mart)
colnames(myannot)[1]="variable"
#add gene lfc
de=read_tsv("DE.genes.tsv")
colnames(de)[1:2]=c("variable","subtype")
topg=merge(top,myannot,by="variable")%>%
	merge(de[,c(1:3,7)],by=c("variable","subtype"))%>%
	arrange(subtype,desc(n))%>%write_tsv("top.tsv")
#get literature on the gene[4] & subtype[2]
known=apply(topg,1,function(x) 
	entrez_search(db="pubmed",
				  term=paste(x[2],x[4],"cancer",sep=' AND ')))
#repeated features
top=top%>%pivot_wider(names_from="subtype",values_from=n)

#which features are most enriched?????????????????????
bp=read_tsv("BP-allFeatures.enrichment")
k=read_tsv("KEGG-allFeatures.enrichment")
enriched=list(bp,k)
names(enriched)=c("BP","KEGG")

features=lapply(enriched,function(x) 
	x%>%dplyr::select(subtype,Description,geneID)%>%
	separate_rows(geneID,sep='/',convert=T))
top=lapply(features,function(x) x%>%group_by(subtype)%>%
	count(geneID)%>%arrange(desc(n))%>%group_map(~head(.x,5)))#top10
top=lapply(top,function(y) do.call(rbind,lapply(1:5,function(x) 
	cbind(subtype=unique(features$BP$subtype)[x],y[[x]]))))
myannot=getBM(attributes = c("hgnc_symbol","entrezgene_id",
	"wikigene_description"),mart=mart)
colnames(myannot)[1]="geneID"
top$BP=merge(top$BP,myannot[,c(1,3)],by="geneID")%>%arrange(subtype)
colnames(myannot)[1:2]=c("hgnc_symbol","geneID")
top$KEGG=merge(top$KEGG,myannot,by="geneID")%>%arrange(subtype)
#didnt add logFC, coz these includes miR & CpGs

#which functions map the same genes for a pair of datasets?????
i=lapply(features,function(x) 
	x%>%distinct(subtype,Description)%>%count(Description)%>%
	filter(n>1)%>%select(Description)%>%unlist)
#use only functions found in at least 2 datsets
features=lapply(1:2,function(x) 
	features[[x]]%>%filter(Description%in%i[[x]]))
#jaccard index
jacc=function(set1,set2){
	inter=length(intersect(set1,set2))
	un=length(union(set1,set2))
return(inter/un)}
intersect_functions=function(data,fun){
	set=data%>%filter(Description==fun)
	subtys=unique(set$subtype)
	#paired contrast
	mat=do.call(rbind,lapply(1:(length(subtys)-1),function(x) 
		do.call(rbind,lapply((1+x):length(subtys),function(y) 
			c("function"=fun,"pair1"=subtys[x],"pair2"=subtys[y],
				"index"=jacc(set$geneID[set$subtype==subtys[x]],
						   set$geneID[set$subtype==subtys[y]]))))))
return(mat)}
jindx=lapply(features,function(x) 
	data.frame(do.call(rbind,lapply(unique(x$Description),function(y)
		intersect_functions(x,y)))))
jindx=lapply(jindx,function(x) cbind(x,"pair"=paste(x$pair1,pair2)))
#plot distributions of Jaccard Index
pdf("enrichJacc.pdf",width=15)
 lapply(jindx,function(z) ggplot(z,aes(x=pair,y=as.numeric(index)))+
 	geom_boxplot()+ylab("jaccard index")+
 	theme(text=element_text(size=18)))
dev.off()
#save it for puma nets
write_tsv(do.call(rbind,lapply(1:2,function(x) 
	cbind("class"=c("BP","KEGG")[x],jindx[[x]][,1:4]))),
"funcJaccI.tsv")


