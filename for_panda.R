###########FEATRUREs MATRIX TO BUILD THE OBJECT NEEDED FOR PANDA
library(tidyverse)
LumB=read_tsv("LumB.selected")
KEGGenrich=read_tsv("KEGG-allFeatures.enrichment")

#components with the function
i=KEGGenrich%>%
  filter(subtype=="LumB"&Description=="MicroRNAs in cancer")%>%
  select(component)%>%unlist
#features in those components
features=unique(LumB$variable[LumB$component%in%i])
#measures of those features
library(data.table)
data=fread("LumB.eigeNormi")
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
data=data[rownames(data)%in%features,]

###########REGULATORY NET TO BUILD THE OBJECT NEEDED FOR PANDA
library(biomaRt)

#regulatory info for CpGs
methy=read_tsv("MapMethy.tsv")
#keep only CpG regulators in feature set
methy=methy[methy$IlmnID%in%features,]
#[1] 1284  of the 1685 selected are llinked to a gene
#keep only interaction with feature set
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","refseq_ncrna", 
	"refseq_mrna"), mart=mart)
myannot=myannot%>%pivot_longer(-1,names_to="type",values_to="refseq")%>%
		filter(refseq!="")
methy=merge(methy,myannot,by="refseq",all.x=T)
methy[methy$ensembl_gene_id%in%features,]

#regulatory info for TFs
tfs=read_tsv("../Downloads/TFtargets.tsv")
tfs=tfs%>%separate_rows(target,sep=',',convert=T)
tfs=tfs%>%separate_rows(TF,sep=',',convert=T)
#keep only interactions between elements of feature set
tfs=tfs[tfs$target%in%features&tfs$TF%in%features,]

#regulatory info for miRNAs
library(multiMiR)
#only hsa-miR-375 has direct targets so u need mature IDs
mirIDs=read_tsv("../Downloads/miR.ids.map.tsv",skip=1)
mirIDs=mirIDs[mirIDs$precursor%in%features,]
#get regulatory interactions with miRNAs in feature set
miRtargets=get_multimir(mirna=c(mirIDs$mature,features),
	summary=F,table="validated",legacy.out=F)
miRtargets=multiMiR::select(miRtargets,keys="validated",
	columns=columns(miRtargets),keytype="type")
#keep only interactions with feature set
miRtargets=miRtargets[miRtargets$target_ensembl%in%features,]
miRtargets=merge(miRtargets,mirIDs,by="mature",all.x=T)

#bring all together
reguEdges=mapply(c, methy[,c("IlmnID","ensembl_gene_id")],
				tfs[,c("TF","target")],
				miRtargets[,c("precursor","target_ensembl")])
colnames(reguEdges)=c("regulator","target")

###########PPI TO BUILD THE OBJECT NEEDED FOR PANDA
library(STRINGdb)
library(igraph)
library(biomaRt)
#get human PPI
string_db <- STRINGdb$new(species=9606)
human_graph <- string_db$get_graph()#Timeout of 60 seconds was reached over & over
#get nodes ensembl_peptide_id
nodes=V(human_graph)$name
nodes=cbind(nodes,gsub("^[0-9]+.","",nodes,perl=T))
colnames(nodes)=c("name","ensembl_peptide_id")
#map ensembl_peptide_id to ensembl_gene_id (in feature set)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=105)
#https://dec2021.archive.ensembl.org
myannot <- getBM(attributes = c("ensembl_gene_id","ensembl_peptide_id"),
     filters = "ensembl_gene_id",values =reguEdges$regulator,
     mart = mart)
#only interested in nodes from PPI that are also TFs in the regulatory net
myannot=merge(myannot,nodes,by="ensembl_peptide_id",all.x=T)
myannot=myannot[!is.na(myannot$name),]
#keep inly the interactions between feature set
edges=get.edgelist(human_graph)
edges=edges[edges[,1]%in%myannot$name&edges[,2]%in%myannot$name,]
#update node names to IDs in feature set 
human_graph=graph.edgelist(edges)
temp=myannot[myannot$name%in%V(human_graph)$name,]
temp=temp[order(match(temp$name,V(human_graph)$name)),]
V(human_graph)$name=temp$ensembl_gene_id
edges=get.edgelist(human_graph)

###########PPI TO BUILD THE OBJECT NEEDED FOR PANDA
#list all together
forPanda=list(expression=data.frame(data),
	motif=data.frame(cbind(reguEdges,1)),
	ppi=data.frame(cbind(edges1,1)))
pandaRes <- panda(forPanda$motif,forPanda$expression,forPanda$ppi,progress=T)
[1] "Initializing and validating"
Error in tfCoopNetwork[Idx] <- ppi[, 3] : 
  NAs are not allowed in subscripted assignments



