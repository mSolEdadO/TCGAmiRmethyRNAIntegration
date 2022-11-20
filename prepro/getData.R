library(SummarizedExperiment)#1.22.0
library(TCGAbiolinks)#2.20.1
library(VennDiagram)#1.6.20

##########SAMPLE IDs PER DATA TYPE#####################
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450")
mthyltn=getResults(mthyltn)
i=substr(mthyltn$cases,1,19)
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type="STAR - Counts")#11/2022 option
#  workflow.type = "HTSeq - Counts")
xprssn=getResults(xprssn)
j=substr(xprssn$cases,1,19)
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification")
mirnas=getResults(mirnas)
k=substr(mirnas$cases,1,19)

##############CONCOURRENT MEASURES########################
sapply(list(i,j,k),function(x) length(unique(x)))
#[1]  890 1217 1202
#[1]  893 1221 1202 #11/2022 
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 848

##############SUBTYPE PER SAMPLE########################
#tissue per sample
samples=data.frame(cbind(samples,
	sapply(samples,function(x) 
		unique(as.character(mthyltn$tissue.definition[i==x])))))
colnames(samples)[2]="tissue"
samples$patient=substr(samples$samples,1,12)

#download subtypes
subtypes=TCGAquery_subtype(tumor="brca")#subtype per patient
sum(!substr(samples,1,12)%in%patients$patient)
#[1] 9
#[1] 3 #11/2022
#only classified samples  are useful
samples=samples[samples$patient%in%subtypes$patient,]
table(subtypes$BRCA_Subtype_PAM50[subtypes$patient%in%samples$patient])
#Basal   Her2   LumA   LumB Normal 
#   129     46    416    140     34 
#normal subtype is discarded NOT normal tissue
samples=samples[which(!samples$patient%in%subtypes$patient[
	subtypes$BRCA_Subtype_PAM50=="Normal"]|
	samples$tissue=="Solid Tissue Normal"),]

#71 patients are duplicated in DB
temp=table(samples[samples$patient%in%samples$patient[duplicated(
	samples$patient)],2:3])
table(apply(temp,2,paste,collapse=""))
#65 patients have 1 tumor & 1 normal samples
#2 patients have 2 tumor samples
#2 patients have 1 tumor & 1 metastatic samples
#2 patients have 1 of each

#subtype per sample
samples$subtype=sapply(samples$patient,function(x) 
	subtypes$BRCA_Subtype_PAM50[subtypes$patient==x])
samples$subtype[samples$tissue=="Solid Tissue Normal"]="Normal"
table(samples$subtype)
# Basal   Her2   LumA   LumB Normal 
#   128     46    416    140     75 
#TCGA-A7-A0CE is missing from Basal coz there's no methylation data
write.table(samples,"subtype.tsv",sep='\t',quote=F,row.names=F)

#plot intersections
i=i[mthyltn$tissue.definition!="Solid Tissue Normal"]
j=j[xprssn$tissue.definition!="Solid Tissue Normal"]
k=k[mirnas$tissue.definition!="Solid Tissue Normal"]
lapply(unique(samples$subtype),function(x)
	venn.diagram(x=list(
		A=unique(j[substr(j,1,12)%in%subtypes$patient[subtypes$BRCA_Subtype_PAM50==x]]),
		B=unique(i[substr(i,1,12)%in%subtypes$patient[subtypes$BRCA_Subtype_PAM50==x]]),
		C=unique(k[substr(k,1,12)%in%subtypes$patient[subtypes$BRCA_Subtype_PAM50==x]])),
		col = "transparent", fill = c("cornflowerblue","green","red"),
		alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
		fontfamily = "sans",cat.fontfamily=rep("sans",3),
		category.names=c("RNAseq","HM450","miRNAseq"),filename=x))

venn.diagram(x=list(
	A=substr(mthyltn$cases[mthyltn$tissue.definition=="Solid Tissue Normal"],1,19),
	B=substr(xprssn$cases[xprssn$tissue.definition=="Solid Tissue Normal"],1,19),
	C=substr(mirnas$cases[mirnas$tissue.definition=="Solid Tissue Normal"],1,19)),
	col = "transparent", fill = c("cornflowerblue","green","red"),
	alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
	fontfamily = "sans",cat.fontfamily=rep("sans",3),
	category.names=c("RNAseq","HM450","miRNAseq"),filename="Normal")

###############ADD CLINICAL INFO##########################
clin <- GDCquery_clinic("TCGA-BRCA","clinical")
clin <- clin[,c("bcr_patient_barcode","gender",
	"tumor_stage","race","vital_status")]
samples=cbind(samples,t(sapply(samples$patient,function(x) 
	clin[clin$bcr_patient_barcode==x,2:4])))
table(subtype$gender)
#female 
#   805 
table(subtype$tumor_stage)
#not reported      stage i     stage ia     stage ib     stage ii    stage iia 
#           6           63           58            4            5          256 
#   stage iib    stage iii   stage iiia   stage iiib   stage iiic     stage iv 
#         196            1          132           20           49           10 
#     stage x 
#           5 
table(subtype$race)
#american indian or alaska native                            asian 
#                               1                               37 
#       black or african american                     not reported 
#                             151                               15 
#                           white 
#                             601 
