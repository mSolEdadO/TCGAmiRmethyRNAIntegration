library(pandaR)

pandaRes <- panda(pandaToyData$motif,#c1:TF,c2:gene with a binding site,c3:1
pandaToyData$expression,#measures
pandaToyData$ppi,#c1:prote 1,c2:prote 2,c3:1
hamming=.1,
progress=TRUE)

###########TO BUILD THE OBJECT NEEDED FOR PANDA U NEED A FEATRUREs MATRIX
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

###########TO BUILD THE OBJECT NEEDED FOR PANDA U NEED A REGULATORY NET
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

##regulatory info for TFs
tfs=read_tsv("../Downloads/TFtargets.tsv")
tfs=tfs%>%separate_rows(target,sep=',',convert=T)
tfs=tfs%>%separate_rows(TF,sep=',',convert=T)
#keep only interactions between elements of feature set
tfs=tfs[tfs$target%in%features&tfs$TF%in%features,]

#CpG to transcript part of motif
library(tidyverse)
library(biomaRt)

methy=fread("../Downloads/MapMethy.tsv")
methy=methy[methy$IlmnID%in%rownames(expre),]
methy=separate_rows(methy[,c("IlmnID","UCSC_RefGene_Name")],
					UCSC_RefGene_Name,convert=T) 
lis