#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))

mi=commandArgs(trailingOnly=TRUE)
edges=read_tsv(mi,col_names=F,show_col_types=F)

g=graph.data.frame(edges[,1:2],directed=F)
E(g)$MI=edges$X3
g=simplify(g,edge.attr.comb="first")#coz MI(a,b)=MI(b,a)
features=V(g)$name

#recover non-duplicated edges
edges=as.data.frame(get.edgelist(g))
edges$MI=E(g)$MI
edges$type=apply(cbind(substr(edges$V1,1,1),substr(edges$V2,1,1)),1,
#so CpG-transcript & transcript-CpG are together
	function(x) paste(sort(x),collapse=''))

#################################real CpG-transcript edges
#library(biomaRt)
methy=read_tsv("/home/msoledad/param-rm//MapMethy.tsv",show_col_types=F)
methy=methy%>%filter(IlmnID%in%features)
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
#	version=105)
#https://dec2021.archive.ensembl.org
#myannot=getBM(attributes = c("ensembl_gene_id","refseq_ncrna", 
#	"refseq_mrna","mirbase_id"), mart=mart)
#write_tsv(myannot,"myannot")
myannot=read_tsv("myannot",show_col_types=F)
methy=myannot%>%pivot_longer(-c(1,4),names_to="type",
	values_to="refseq")%>%merge(methy,by="refseq",all.y=T)
methy=methy%>%select(IlmnID,ensembl_gene_id,mirbase_id)%>%
	pivot_longer(-1,values_to="target",names_to="type")%>%
	filter(!is.na(target))
#get MI between every pair. 
#how does it compare with obtained CpG-transcript edges?

#################################real TF-transcript edges
tfs=read_tsv("/home/msoledad/param-rm/TFtargets.tsv",show_col_types=F)
tfs=tfs%>%separate_rows(target,sep=',',convert=T)%>%
	separate_rows(TF,sep=',',convert=T)%>%
	filter(target%in%features|TF%in%features)

#################################real miR-transcript edges
if(length(grep("hsa",features))>0){
suppressPackageStartupMessages(library(multiMiR))
mirIDs=read_tsv("/home/msoledad/param-rm/miR.ids.map.tsv",skip=1,show_col_types=F)
#add precursors as possible regulators
temp=mirIDs[,c(1,1)]
colnames(temp)[2]="mature"
temp$mature=gsub("mir","miR",temp$mature)
mirIDs=unique(rbind(mirIDs,temp))
mirIDs=mirIDs[mirIDs$precursor%in%features,]
#get regulatory interactions with miRNAs in feature set
miRtargets=get_multimir(mirna=c(mirIDs$mature,features),
	summary=F,table="validated",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="validated",
	columns=columns(miRtargets),keytype="type")
colnames(miRtargets)[3]="mature"
miRtargets=merge(miRtargets,mirIDs,by="mature")
#join all pairs
reguEdges=mapply(c, methy[,c("IlmnID","target")],
				tfs[,c("TF","target")],
				miRtargets[,c("precursor","target_ensembl")])
}else{
print("No miRNAs")	
reguEdges=mapply(c, methy[,c("IlmnID","target")],
				tfs[,c("TF","target")])
}
gr=graph.data.frame(reguEdges)

if(nrow(reguEdges)>0){
###########################################get real MI values
library(infotheo)

subtype=unlist(strsplit(mi,".",fixed=T))[2]
#expression/methylation data
data=suppressWarnings(data.table::fread(paste(subtype,"eigeNormi",sep='.')))
data=data[data$V1%in%V(gr)$name,]
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)

#keep just the interactions involving nodes in the data
gr=suppressWarnings(subgraph(gr,which(V(gr)$name%in%rownames(data))))
reguEdges=get.edgelist(gr)

#methods used neeed discrete data
discrData=discretize(t(data))
#can chage it for parLapply if reguEdges is too large
reguEdges=data.frame(cbind(reguEdges,apply(reguEdges,1,function(x) 
	condinformation(discrData[,x[1]],discrData[,x[2]]))))
colnames(reguEdges)=c("regulator","target","MI")
reguEdges$type=apply(cbind(substr(reguEdges$regulator,1,1),
						substr(reguEdges$target,1,1)),1,
#so CpG-transcript & transcript-CpG are together
				function(x) paste(sort(x),collapse=''))
#####################################compare with MI in g
#keep only MI values>median(MI for real regulatory interactions)
thres=sapply(unique(reguEdges$type),function(x) 
	median(as.numeric(reguEdges$MI[reguEdges$type==x])))
filtered=do.call(rbind,lapply(names(thres),function(x) 
				edges%>%filter(type==x&MI>=thres[x])))
write_tsv(filtered,file=gsub("sort","filtered",mi))
}
#just gonna ignore CpG-CpG pair?????????????????YEP
edges1=edges%>%filter(type!="cc")

#is the MI distri the same across omics?
#ggplot(edges1,aes(x=MI,col=type))+geom_density(aes(y=..scaled..))+scale_x_continuous(trans="log10")
pvalMat=sapply(unique(edges1$type),function(x) 
	sapply(unique(edges1$type),function(y) 
		suppressWarnings(ks.test(edges1$MI[edges1$type==x],
			edges1$MI[edges1$type==y]))$p.val))
diffTypes=which(pvalMat<0.05,arr.ind=T)
if(length(diffTypes)>0){
	i=nrow(diffTypes)/2
	sapply(i,function(x) paste(colnames(pvalMat)[diffTypes[x,1]],
		"is different than",colnames(pvalMat)[diffTypes[x,2]]))
}else{
#if MI distributions were the same, u could use DPI
#or choose the threshold from any type
	filtered=edges1%>%filter(MI>=min(thres))
	write_tsv(filtered,file=gsub("sort","filtered.alt",mi))
}
temp=as.data.frame(mapply(c,cbind("predi",edges1),
	cbind("known",reguEdges),
	cbind("filtered",filtered)))
table(temp[,c(1,5)])
