library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)

#get sample ids per data type
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450")
mthyltn=getResults(mthyltn)
i=substr(mthyltn$cases,1,15)
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts")
xprssn=getResults(xprssn)
j=substr(xprssn$cases,1,15)
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification")
mirnas=getResults(mirnas)
k=substr(mirnas$cases,1,15)

sapply(list(i,j,k),function(x) length(unique(x)))
#[1]  890 1217 1202
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 848

#download subtypes
subtype=TCGA_MolecularSubtype(xprssn$cases)$subtypes
#https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018
subtype$barcode=substr(subtype$samples,1,15)
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

#some patients were measured several times
sum(duplicated(subtype$patients))
#[1] 71
#12 are NOT normal-tumor pairs
temp=table(subtype[subtype$patients%in%
	subtype$patients[duplicated(subtype$patients)],2:3])
sum(temp[5,]!=1)
#[1] 12
#TCGA-E9-A1NF has 1 LumA sample & 1 LumB
#TCGA-BH-A1FE has 2 normal samples + 1 LumA

#can duplicates be explained by metastasis?
#sum(subtype$samples%in%xprssn$cases)
tissue=xprssn[xprssn$cases%in%subtype$samples,c(4,28)]
tissue=tissue[order(match(tissue$cases,subtype$samples)),]
subtype=cbind(subtype,tissue$tissue.definition)
colnames(subtype)[5]="tissue"
#problems with sample IDs
lapply(unique(subtype$tissue),function(x) 
	table(subtype[subtype$tissue==x,c(2,5)]))