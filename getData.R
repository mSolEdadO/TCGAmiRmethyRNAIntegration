library(TCGAbiolinks)
library(VennDiagram)

#indexes of each data type
xprssn <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	sample.type = "Primary solid Tumor",
	workflow.type = "HTSeq - Counts")
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform = "Illumina Human Methylation 450",
	sample.type = "Primary solid Tumor")
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor")

#unique patients for each data
patients<-function(indexes){
	substr(getResults(indexes)$cases,1,12)}
xC=unique(patients(xprssn))
#length(xC)
#[1] 1091 
miC=unique(patients(mirnas))
#length(miC)
#[1] 1078
mtC=unique(patients(mthyltn))
#length(mtC)
#[1] 1095

#patients intersecting the three data types
venn.diagram(x = list(A=xC,B=miC,C=mtC), filename = "Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green","red"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="white", 
cat.cex = 1.5,margin = 0.1,category.names=c("expression","miRNAs",
"methylation"))
inter=intersect(intersect(xC,mtC),miC)

#data of each type from patients in the intersection
xC=getResults(xprssn)$cases[patients(xprssn)%in%inter])
xprssn <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	sample.type = c("Primary solid Tumor"),
	workflow.type = "HTSeq - Counts",barcode=xC)
#GDCdownload(xprssn) 
miC=getResults(mirnas)$cases[patients(mirnas)%in%inter]
mirna <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor",barcode=miC)
#GDCdownload(mirna)
mtC=getResults(mthyltn)$cases[patients(mthyltn)%in%inter]
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	sample.type = c("Primary solid Tumor"),barcode=mtC)
#GDCdownload(mthyltn)

#Data from normal tissue
xprssN <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	sample.type = c("Solid Tissue Normal"),
	workflow.type = "HTSeq - Counts")
mthyltN <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform = "Illumina Human Methylation 450",
	sample.type = c("Solid Tissue Normal"))
mirnasN <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Solid Tissue Normal")

#Tumor and normal tissue matched samples
venn.diagram(x = list(A=xC,B=unique(patients(xprssN))), filename = "ExpressionTN.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="white", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))