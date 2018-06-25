#find patiets tested for all the platforms for tumor and normal tissue
interX=intersect(unique(patients(xprssn)),unique(patients(xprssN)))
interMi=intersect(unique(patients(mirna)),unique(patients(mirnasN)))
interMT=intersect(unique(patients(mthyltn)),unique(patients(mthyltN)))
inter=intersect(intersect(interX,interMi),interMT)

#get sample ids for each platform
xC=unique(c(getResults(xprssn)$cases[patients(xprssn)%in%inter],
			getResults(xprssN)$cases[patients(xprssN)%in%inter]))
miC=unique(c(getResults(mirna)$cases[patients(mirna)%in%inter],
			getResults(mirnasN)$cases[patients(mirnasN)%in%inter]))
mtC=unique(c(getResults(mthyltn)$cases[patients(mthyltn)%in%inter],
			getResults(mthyltN)$cases[patients(mthyltN)%in%inter]))

#get the data
xprssn <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "Gene Expression Quantification",
	workflow.type = "HTSeq - Counts",
	barcode=xC)
GDCdownload(xprssn)
#GDCdownload will download 190 files. A total of 48.469537 MB
mirna <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	barcode=miC)
GDCdownload(mirna)
#GDCdownload will download 191 files. A total of 9.594273 MB
mthyltn <- GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	barcode=mtC)
GDCdownload(mthyltn)
#GDCdownload will download 243 files. A total of 26.450998939 GB
