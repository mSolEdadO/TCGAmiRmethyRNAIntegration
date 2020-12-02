library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)

#get sample ids per data type
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450")
j=substr(getResults(mthyltn)$cases,1,12)
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts")
i=substr(getResults(xprssn)$cases,1,12)
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification")
k=substr(getResults(mirnas)$cases,1,12)

sapply(list(j,i,k),function(x) length(unique(x)))
#[1]  789 1092 1079
#matching ids
patients=intersect(j,intersect(i,k))
length(patients)
#[1] 777

#sample subtypes
subtype=TCGA_MolecularSubtype(
	getResults(xprssn)$cases)$subtypes
sum(!patients%in%subtype$patients)#non classified samples
#[1] 5
#plot intersections
lapply(unique(subtype$subtype),function(x)
	venn.diagram(x=list(A=unique(i[i%in%subtype$patients[subtype$subtype==x]]),
		B=unique(j[j%in%subtype$patients[subtype$subtype==x]]),
		C=unique(k[k%in%subtype$patients[subtype$subtype==x]])),
		col = "transparent", fill = c("cornflowerblue","green","red"),
		alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
		fontfamily = "sans",cat.fontfamily=rep("sans",3),
		category.names=c("RNAseq","HM450","miRNAseq"),filename=x))

#get data of classified samples only 
patients=patients[patients%in%subtype$patients]
table(unique(subtype[,c(2,4)])[,1])#do not count duplications
# BRCA.Basal   BRCA.Her2   BRCA.LumA   BRCA.LumB BRCA.Normal 
#        127          46         420         146         112 

#get ALL sample ids for these patients
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450",
	barcode=patients)
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  barcode=patients)
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	barcode=patients)