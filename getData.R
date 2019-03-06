library(TCGAbiolinks)
library(VennDiagram)

#indexes of each data type
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	sample.type = "Primary solid Tumor",
	platform="Illumina Human Methylation 450")
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor")
#unique patients for each data
miC=unique(substr(getResults(mirnas)$cases,1,12))
length(miC)
#[1] 1013
mtC=unique(substr(getResults(mthyltn)$cases[getResults(mthyltn)$platform=="Illumina Human Methylation 450"],1,12))
length(mtC)
#[1] 782

#intersection with the other data types
venn.diagram(x = list(A=unique(substr(subtipos$patient[subtipos$pbcmc2=="Her2"],1,12)),B=miC,C=mtC),
 filename = "Her2-450.tiff",col = "transparent", fill = c("cornflowerblue","green","red"),
 alpha = 0.50,cex = 1.5, fontface = "bold", label.col="white", cat.cex = 1.5,margin = 0.1,
 category.names=c("expression","miRNAs","methylation"),main="Her2-450")
venn.diagram(x = list(A=unique(substr(subtipos$patient[subtipos$pbcmc2=="Her2"],1,12)),B=miC,C=unique(patients(mthyltn))),
 filename = "Her2.tiff",col = "transparent", fill = c("cornflowerblue","green","red"), alpha = 0.50,
 cex = 1.5, fontface = "bold", label.col="white", cat.cex = 1.5,margin = 0.1,
 category.names=c("expression","miRNAs","methylation"),main="Her2")
venn.diagram(x = list(A=unique(substr(subtipos$patient[subtipos$pbcmc2=="LumA"],1,12)),B=miC,C=mtC),
 filename = "LumA.tiff",col = "transparent", fill = c("cornflowerblue","green","red"),
  alpha = 0.50,cex = 1.5, fontface = "bold", label.col="white", cat.cex = 1.5,margin = 0.1,
  category.names=c("expression","miRNAs","methylation"),main="LumA")
venn.diagram(x = list(A=unique(substr(subtipos$patient[subtipos$pbcmc2=="LumB"],1,12)),B=miC,C=mtC), 
  filename = "LumB.tiff",col = "transparent", fill = c("cornflowerblue","green","red"), alpha = 0.50,
  cex = 1.5, fontface = "bold", label.col="white", cat.cex = 1.5,margin = 0.1,
  category.names=c("expression","miRNAs","methylation"),main="LumB")
venn.diagram(x = list(A=unique(substr(subtipos$patient[subtipos$pbcmc2=="Basal"],1,12)),B=miC,C=mtC), 
  filename = "Basal.tiff",col = "transparent", fill = c("cornflowerblue","green","red"), alpha = 0.50,
  cex = 1.5, fontface = "bold", label.col="white", cat.cex = 1.5,margin = 0.1,
  category.names=c("expression","miRNAs","methylation"),main="Basal")
#convert -append Basal.tiff LumB.tiff LumA.tiff Her2-450.tiff Her2.tiff subtipos.tiff
intersec=sapply(c("LumA","LumB","Basal","Her2"),function(x) 
	intersect(intersect(unique(substr(subtipos$patient[subtipos$pbcmc2==x],1,12)),mtC),miC))
#intersec$Her2=intersect(intersect(unique(substr(subtipos$patient[subtipos$pbcmc2=="Her2"],1,12)),
#	unique(substr(getResults(mthyltn)$cases,1,12))),miC)
#sapply(intersec,length)
# LumA  LumB Basal  Her2 
#  331   178   135   119 

#data of each type from patients in the intersection
miC=getResults(mirnas)$cases[substr(getResults(mirnas)$cases,1,12)%in%unlist(intersec)]
#length(miC)
#[1] 733
mirna <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Primary solid Tumor",barcode=miC)
GDCdownload(mirna)
mtC=getResults(mthyltn)$cases[substr(getResults(mthyltn)$cases,1,12)%in%unlist(intersec)]
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	sample.type = c("Primary solid Tumor"),barcode=mtC)
GDCdownload(mthyltn)

#Data from normal tissue
mthyltN <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform = "Illumina Human Methylation 450",
	sample.type = c("Solid Tissue Normal"))
mirnasN <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	sample.type = "Solid Tissue Normal")
