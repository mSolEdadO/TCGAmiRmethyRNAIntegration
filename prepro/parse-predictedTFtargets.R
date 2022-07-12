library(tftargets)
library(biomaRt)
library(tidyverse)

#############################TABLE FROM LISTS
#entrez ID for targets
#TRED: Predicted and known human transcription factor targets.
TRED=as.data.frame(do.call(rbind,lapply(1:length(TRED),
	function(x) cbind(TF=names(TRED)[x],target=TRED[[x]]))))

#gene symbols
#ITFP: Predicted human transcription factor targets 
ITFP=as.data.frame(do.call(rbind,lapply(1:length(ITFP),
	function(x) cbind(TF=names(ITFP)[x],target=ITFP[[x]]))))

#per tissue w gene symbols
#Marbach2016: close to 400 cell type- and tissue-specific gene regulatory networks for human
Marbach2016=as.data.frame(do.call(rbind,lapply(1:length(Marbach2016),
	function(x) cbind(TF=names(Marbach2016)[x],target=Marbach2016[[x]]))))
#since columns are same
predicted=rbind(cbind(db="ITFP",ITFP),
	cbind(db="Marbach",Marbach2016))

############################GENE IDS
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"hgnc_symbol"),mart=mart)

#all IDs for entrez targets
colnames(myannot)[2]="target"
TRED=merge(TRED,myannot,by="target",all.x=T)
colnames(TRED)[1]=c("entrez")
TRED$db="TRED"

#all IDs for symbol targets
colnames(myannot)[2:3]=c("entrez","target")
predicted=merge(predicted,myannot,by="target",all.x=T)
colnames(predicted)[1]="hgnc_symbol"

#paste all together
predicted=rbind(TRED[,colnames(predicted)],predicted)
colnames(predicted)[4]="target_ensembl"

#drop targets without ensembl, coz TCGA has ensembl
predicted=predicted%>%filter(!is.na(ensembl_gene_id))
#get TF ensembl
colnames(predicted)[4]="target_ensembl"
colnames(myannot)[3]="TF"
predicted=merge(predicted,myannot[,c(1,3)],by="TF",all.x=T)
colnames(predicted)[6]="tf_ensembl"
predicted=predicted%>%filter(!is.na(tf_ensembl))
#I've made a monster
#dim(predicted)
#[1] 1398888       6
write_tsv(tfs,file="prediTFtargets.tsv")
