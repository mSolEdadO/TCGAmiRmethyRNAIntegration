#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
mtrx=commandArgs(trailingOnly=TRUE)

library(biomaRt)
library(tidyverse)
library(multiMiR)
library(STRINGdb)
library(igraph)

###########FEATRURE MATRIX TO BUILD THE OBJECT NEEDED BY PUMA
data=data.table::fread(mtrx)
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
features=rownames(data)

###########METHY HAS NO ALTERNATIVE

###########PREDICTED TF-TARGETS
tfs=read_tsv("/home/msoledad/param-rm/prediTFtargets.tsv")
#keep only TFs in the feature set
tfs=tfs[tfs$tf_ensembl%in%features,]
#keep only targets in the feature set
i=tfs[tfs$target_ensembl%in%features,c("tf_ensembl","target_ensembl")]

#recover TFs targeting miRNAs too
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
#	version=105)
#https://dec2021.archive.ensembl.org
#myannot=getBM(attributes = c("ensembl_gene_id","mirbase_id"), mart=mart)
#write_tsv(myannot,"myannot")
myannot=read_tsv("myannot")
myannot=myannot[myannot$mirbase_id%in%features,]
j=tfs[tfs$target_ensembl%in%myannot$ensembl_gene_id,
      c("tf_ensembl","target_ensembl")]
colnames(j)[2]="ensembl_gene_id"
j=merge(j,myannot,by="ensembl_gene_id")
tfs=data.frame(mapply(c,i,j[,c("tf_ensembl","mirbase_id")]))
ft=tryCatch(nrow(tfs),error=0)

##############miRNAs
#only hsa-miR-375 has direct targets so u need mature IDs
mirIDs=read_tsv("/home/msoledad/param-rm/miR.ids.map.tsv",skip=1)
#add precursors as possible regulators
temp=mirIDs[,c(1,1)]
colnames(temp)[2]="mature"
temp$mature=gsub("mir","miR",temp$mature)
mirIDsm=unique(rbind(mirIDs,temp))
mirIDs=mirIDs[mirIDs$precursor%in%features,]
#get regulatory interactions with miRNAs in feature set
miRtargets=get_multimir(mirna=c(mirIDs$mature,features),
	summary=F,table="predicted",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="predicted",
	columns=columns(miRtargets),keytype="type")
#keep only interactions with feature set
miRtargets=miRtargets[miRtargets$target_ensembl%in%features,]
fmi=tryCatch(nrow(miRtargets),error=0)

miRtargets=if(fmi>0){
	colnames(miRtargets)[3]="mature"
	merge(miRtargets,mirIDs,by="mature",all.x=T)}

#bring all together
if(fmi==0 & ft==0){stop("No prior")}

reguEdges=mapply(c,tfs,
			    miRtargets[,c("precursor","target_ensembl")])
#reguEdges=reguEdges[rowSums(apply(reguEdges,2,function(x) x=="NA"))==0,]#drop NAs
if(class(reguEdges)=="character"){reguEdges=t(data.frame(reguEdges))}
write_tsv(as.data.frame(cbind(reguEdges,1)),
	gsub("mtrx","predi.moti",mtrx),col_names=F)
######################NON TF LIST NEEDED BY PUMA
writeLines(features[!features%in%tfs$TF],
	gsub("mtrx","nonTF",mtrx))

###########PPI TO BUILD THE OBJECT NEEDED BY PUMA
if(ft<1){stop("No PPI")}
#get human PPI
#string_db <- STRINGdb$new(species=9606,score_threshold=999)
#human_graph <- string_db$get_graph()#Timeout of 60 seconds was reached over & over
#human_edges=get.edgelist(human_graph)
#write_tsv(as.data.frame(human_edges),"human_graph")
#loading pre-saved net avoids crashes due to Timeout but don't get updated!!!!
human_edges=read_tsv("human_graph")
human_graph=graph.data.frame(human_edges,directed=F)

#get nodes ensembl_peptide_id
nodes=V(human_graph)$name
nodes=cbind(nodes,gsub("^[0-9]+.","",nodes,perl=T))
colnames(nodes)=c("name","ensembl_peptide_id")
#map ensembl_peptide_id to ensembl_gene_id (in feature set)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=105)
myannot <- getBM(attributes = c("ensembl_gene_id","ensembl_peptide_id"),
     mart = mart)
write_tsv(myannot,"myannot.Protes")
#myannot=read_tsv("myannot.Protes")                                         
myannot=myannot[myannot$ensembl_gene_id%in%tfs$TF,]
#only interested in nodes from PPI that are also TFs in the regulatory net
myannot=merge(myannot,nodes,by="ensembl_peptide_id",all.x=T)
myannot=myannot[!is.na(myannot$name),]
#keep only the interactions between feature set
human_edges=human_edges[human_edges$V1%in%myannot$name&
	human_edges$V2%in%myannot$name,]
#update node names to IDs in feature set 
human_graph=graph.data.frame(human_edges)
temp=myannot[myannot$name%in%V(human_graph)$name,]
temp=temp[order(match(temp$name,V(human_graph)$name)),]
V(human_graph)$name=temp$ensembl_gene_id
edges=get.edgelist(human_graph)
#edges=edges[rowSums(apply(edges,2,function(x) x=="NA"))==0,]#drop NAs
if(nrow(edges)>0){
	write_tsv(cbind(as.data.frame(edges),1),
		gsub("mtrx","pp",mtrx),col_names=F)}

print(paste(nrow(reguEdges),"regulatory interactions",
						nrow(edges),"ppi between TFs",sep=' '))

#python run_puma.py -e expre -m moti -p pp -i nonTF.txt -o mipuma.txt
##############################
#library(igraph)
#library(tidyverse)
##library(ggplot2)
##library(RCy3)

#puma=read_tsv("temp",col_names=F)
#temp=cbind(substr(puma$X1,1,1),substr(puma$X2,1,1))
#puma$type=apply(temp,1,function(x) paste(sort(x),collapse='-'))
#g=graph.data.frame(puma[,1:2])#puma produces DIRECTED graphs
#E(g)$weight=puma$X4
#E(g)$type=puma$type
##png("cor.png")
##ggplot(puma,aes(x=abs(X4),fill=type,color=type))+
##geom_density(aes(y=..scaled..),alpha=0.4)+
##scale_x_continuous(trans="log10")
##dev.off()
#g1=subgraph.edges(g,which(abs(E(g)$weight)>1))
#V(g1)$type=substr(V(g)$name,1,1)
#V(g1)$size=degree(g1)
#gc=createNetworkFromIgraph(g1,"myIgraph")

#mi=read_tsv("LumB.miRNAsCancer.MI")
#temp=cbind(substr(mi$V1,1,1),substr(mi$V2,1,1))
#mi$type=apply(temp,1,function(x) paste(sort(x),collapse='-'))
#g=graph.data.frame(mi[,1:2],directed=F)#puma produces DIRECTED graphs
#known=read_tsv("moti",col_names=F)
#known$type=paste(substr(known$X1,1,1),substr(known$X2,1,1),sep="-")
##c-E E-E h-E 
## 26  10 564 
#A=g[known$X1,known$X2]
#known$X3=sapply(1:nrow(A),function(x) A[x,x])
#ks.test(mi$MI[mi$type=="c-E"],known$X3[known$type=="c-E"])
##D = 0.47164, p-value = 1.899e-05
#ks.test(mi$MI[mi$type=="E-E"],known$X3[known$type=="E-E"])
##D = 0.3, p-value = 0.3308
#ks.test(mi$MI[mi$type=="E-h"],known$X3[known$type=="h-E"])
	#D = 0.058905, p-value = 0.07386
