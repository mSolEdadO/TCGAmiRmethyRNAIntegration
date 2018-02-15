#FALTA VER COMO BAJAR EL ARCHIVO CORECTO PARA EL SUMARY:
# "results"
library(TCGAbiolinks)

#get Gene Expresion data
query <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	sample.type = c("Primary solid Tumor"))
GDCdownload(query) 

#get miRNA Expresion data
query <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor")
GDCdownload(query)

#get methylation data
query <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	sample.type = c("Primary solid Tumor"),
	legacy = TRUE)
GDCdownload(query)
# con barcode puedes bajar tar.gz
