library(TCGAbiolinks)
library(SummarizedExperiment)
library(NOISeq)
library(Venn.Diagram)
library(doParallel)

#####################################prepare the data########################################
#reference data
mthyltN=GDCprepare(mthyltN)
temp=cbind(colnames(mthyltN),"Illumina Human Methylation 450",substr(colnames(mthyltN),1,12),"normal")
methyDesign=rbind(methyDesign,temp)
mthyltN=unlist(assays(mthyltN))
mthyltN=mthyltN[nas<ncol(mthyltN)*0.75,]#me quedo solo con las probes con dato para al menos 75% de las muestras
#tumor 450 platform data
save(methyDesign,mthyltN,file="GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy_ori.RData")
mtC=methyDesign[methyDesign[,4]!="Her2"&methyDesign[,4]!="normal",1]
mthyltn=GDCquery(project = "TCGA-BRCA",data.category = "DNA Methylation",barcode=mtC)
methy450=GDCprepare(mthyltn)
annotMethy=rowData(methy450)
methy450=unlist(assays(methy450))
nas=rowSums(is.na(methy450))
#table(nas) hay 365411 probes completas y 89512 sin dato
methy450=methy450[nas<ncol(methy450),]
save(annotMethy,methy450,file="GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy450.RData")
#tumor 27 platform data
#mtC=methyDesign[methyDesign[,4]=="Her2",1]
#mthyltn=GDCquery(project = "TCGA-BRCA",data.category = "DNA Methylation",barcode=mtC)
#methy27=GDCprepare(mthyltn)
#dim: 25978 119 
#methy27=unlist(assays(methy27))
#nas=rowSums(is.na(methy27))#hay 22746 probes completas y 119 sin dato
#methy27=methy27[nas<ncol(methy27),]
#save(methy27,file="GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy27.RData")
##############################impute missing data########################################
load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy450.RData")
load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methyN_ori.RData")
#load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy27.RData")
methySubti=lapply(c("Basal","LumA","LumB","Her2"),function(x) 
	methy450[,colnames(methy450)%in%methyDesign[methyDesign[,4]==x,1]])
methySubti$normal=mthyltN
#methySubti$Her2=methy27
sapply(methySubti,dim)
#      Basal   LumA   LumB normal   Her2
#[1,] 396065 396065 396065 395814 396065
#[2,]    137    337    178     96     75
nas=lapply(methySubti,function(x) rowSums(is.na(x)))
methySubti=lapply(1:5,function(x) methySubti[[x]][nas[[x]]<ncol(methySubti[[x]])*0.75,])
sapply(methySubti,dim)
#      Basal   LumA   LumB normal   Her2
#[1,] 395820 395818 395817 395814 395815
#[2,]    137    337    178     96    75
#impute missing values
no_cores=detectCores()-1
registerDoParallel(cores=no_cores)
cl=makeCluster(no_cores,type="FORK")
methySubti=parLapply(cl,methySubti,function(x) impute.knn(x,k=15,maxp=nrow(x))$data)
names(methySubti)=c("Basal","LumA","LumB","Her2","normal")
#Buuren and Groothuis suggested that 15 variables are generally sufficient for imputation
#maxp forces to impute always with knn, even when there's too much NA in columns/rows
stopCluster(cl)
save(methySubti,methyDesign,file="imptMethy.RData")

#check the intersection between normal & tumor samples
venn.diagram(x = list(A=methyDesign$patient[methyDesign$shortLetterCode!="NT"],
	B=methyDesign$patient[methyDesign$shortLetterCode=="NT"]), filename = "matchMethyNormalTumor.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))
####################################correct batches#########################################
colnames(methyDesign)=c("barcode","platform","patient","subtype")
methyDesign=as.data.frame(methyDesign)

probes=intersect(intersect(intersect(intersect(rownames(methySubti$Basal),rownames(methySubti$LumA)),rownames(methySubti$LumB)),rownames(methySubti$Her2)),rownames(methySubti$normal))
methySubti=lapply(methySubti,function(x) x[rownames(x)%in%probes,])
#Beta values have severe heteroscedasticity for highly methylated or unmethylated CpG sites.
#M-value method provides much better performance in terms of detection rate and true positive rate
#for both highly methylated and unmethylated CpG sites
Msubti=lapply(methySubti,function(x) log2(x/(1-x)))
noiseqData=lapply(Msubti,function(x) readData(data = x, factors = methyDesign[methyDesign$barcode%in%colnames(x),]))
myPCAsubti = lapply(noiseqData,function(x) dat(x, type = "PCA", norm = T,logtransf=T))
shared=do.call(cbind,Msubti)
noiseqData1=readData(data = shared, factors = methyDesign[methyDesign$barcode%in%colnames(shared),])
myPCAshared=dat(noiseqData1, type = "PCA", norm = T,logtransf=T)
pdf("mPCA.pdf")
par(mfrow=c(3,2) )
explo.plot(myPCAsubti$Basal, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCAsubti$LumA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCAsubti$LumB, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCAsubti$Her2, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCAsubti$normal, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCAshared, samples = c(1,2), plottype = "scores", factor = "subtype")
dev.off()

methySubti=lapply(methySubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])#keep only one sample per patient
#sapply(methySubti,dim)
#       LumA  Basal   LumB normal
#[1,] 395808 395808 395808 395808
#[2,]    331    135    177     96
save(methySubti,methyDesign,Msubti,file="imptMethy.RData")
