library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)

#get sample ids per data type
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450")
mthyltn=getResults(mthyltn)
i=substr(mthyltn$cases,1,19)
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts")
xprssn=getResults(xprssn)
j=substr(xprssn$cases,1,19)
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification")
mirnas=getResults(mirnas)
k=substr(mirnas$cases,1,19)

sapply(list(i,j,k),function(x) length(unique(x)))
#[1]  890 1217 1202
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 848

#sample subtypes
subtype=TCGA_MolecularSubtype(xprssn$cases)$subtypes
subtype$barcode=substr(subtype$samples,1,19)
sum(!samples%in%subtype$barcode)#non classified samples
#[1] 6
#only classified samples  are useful
samples=samples[samples%in%subtype$barcode]
table(unique(subtype[subtype$barcode%in%samples,c(1,2)])[,2])#do not count duplicates
# BRCA.Basal   BRCA.Her2   BRCA.LumA   BRCA.LumB BRCA.Normal 
#        126          46         423         146         101 

#plot intersections
lapply(unique(subtype$subtype),function(x)
	venn.diagram(x=list(
		A=unique(j[j%in%subtype$barcode[subtype$subtype==x]]),
		B=unique(i[i%in%subtype$barcode[subtype$subtype==x]]),
		C=unique(k[k%in%subtype$barcode[subtype$subtype==x]])),
		col = "transparent", fill = c("cornflowerblue","green","red"),
		alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
		fontfamily = "sans",cat.fontfamily=rep("sans",3),
		category.names=c("RNAseq","HM450","miRNAseq"),filename=x))

subtype=subtype[subtype$barcode%in%samples,]
write.table(subtype,"subtype.tsv",sep='\t',quote=F,row.names=F)

#some samples were measured twice
sum(duplicated(i[i%in%samples]))
#[1] 2
sum(duplicated(j[j%in%samples]))
#[1] 4
sum(duplicated(k[k%in%samples]))
#[1] 4
