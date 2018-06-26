library(genefu)
library(biomaRt)  
library(TCGAbiolinks)

#load matrixes done with checkMatched-mRNA.R
tumor=read.table("GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/mproteTumor.mtrx")
normal=read.table("GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/mproteNormal.mtrx")

#get entrez IDs for genefu
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "entrezgene"),values=rownames(tumor),
 mart=mart)
colnames(myannot)=c("probe","EntrezGene.ID")#genefu needs this colnames

#make IDs list and matrix concide
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
data=cbind(tumor,normal)
data=data[rownames(data)%in%myannot$ensembl_gene_id,]
data=t(data)#genefu needs probes in colums and samples in rows
 
#subtyping
subtypes=molecular.subtyping(sbt.model="pam50",data=data,annot=myannot,do.mapping=T)
#table(subtypes$subtype)
# Basal   Her2   LumB   LumA Normal 
#   147     87    321    217    109 

#confussion matrix
lala=GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",workflow.type = "HTSeq - Counts",
	barcode=rownames(data))
lala=GDCprepare(lala,directory="GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/")
knownSubtypes=colData(lala)
rm(lala)
confussion=knownSubtypes[,c(2,77)]
confussion=cbind(as.data.frame(confussion),unlist(subtypes$subtype))
table(confussion[,2:3])
                  predicted
subtype_PAM50.mRNA Basal Her2 LumB LumA Normal <NA>
     Basal-like        4    4   14    9      4    0
     HER2-enriched     3    0    4    3      4    0
     Luminal A        13   11   51   20     11    0
     Luminal B         7    5   13   15      6    0
     Normal-like       0    0    2    3      0    0
     <NA>            120   67  237  167     84    0