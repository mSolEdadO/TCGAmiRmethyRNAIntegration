library(TCGAbiolinks)

#get gene expresion indexes
xprssn <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	sample.type = c("Primary solid Tumor"),
	workflow.type = "HTSeq - Counts")#archivo correcto?
#get methylation indexes
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	sample.type = c("Primary solid Tumor"))
#get miRNA Expresion indexes
mirna <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor")

#unique patients with each data
patients<-function(indexes){
	substr(getResults(indexes)$cases,1,12)}
xC=unique(patients(xprssn))
#length(xC)
#[1] 1091 
miC=unique(patients(mirnas))
#length(miC)
#[1] 1078
mtC=unique(patients(mthyltn))
#length()
#[1] 1095
a=venn(list(xC,mtC,miC),show.plot=FALSE)
#xprssn=getResults(xprssn)[patients(xprssn)%in%attr(a,"intersections")$'A:B:C',])
#mirnas=getResults(mirnas)[patients(mirnas)%in%attr(a,"intersections")$'A:B:C',])
#mthyltn=getResults(mthyltn)[patients(mthyltn)%in%attr(a,"intersections")$'A:B:C',])


#dowload data of patients with the three studies
GDCdownload(xprssn) 
GDCdownload(mirnas)
GDCdownload(mthyltn)