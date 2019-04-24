library(TCGAbiolinks)
library(SummarizedExperiment)
library(Venn.Diagram)
library(doParallel)

#######prepare the data########################################
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
#annotMethy=rowData(methy450)
methy450=unlist(assays(methy450))
nas=rowSums(is.na(methy450))
#table(nas) hay 365411 probes completas y 89512 sin dato
methy450=methy450[nas<ncol(methy450),]
save(annotMethy,methy450,file="GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy450.RData")
save(methy450,file="GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy450.RData")
#######impute missing data########################################
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
############filter noise-prone probes######################################
library(lumi)
library(minfi)

probes=intersect(intersect(intersect(intersect(rownames(methySubti$Basal),rownames(methySubti$LumA)),rownames(methySubti$LumB)),rownames(methySubti$Her2)),rownames(methySubti$normal))
methySubti=lapply(methySubti,function(x) x[rownames(x)%in%probes,])

methyDesign=as.data.frame(methyDesign)
methyDesign=cbind(methyDesign,sapply(strsplit(as.character(methyDesign$barcode),"-"),function(x) substr(x[4],1,2)))
methyDesign=cbind(methyDesign,sapply(strsplit(as.character(methyDesign$barcode),"-"),function(x) substr(x[4],3,3)))
methyDesign=cbind(methyDesign,sapply(strsplit(as.character(methyDesign$barcode),"-"),function(x) substr(x[5],1,2)))
methyDesign=cbind(methyDesign,sapply(strsplit(as.character(methyDesign$barcode),"-"),function(x) substr(x[5],3,3)))
methyDesign=cbind(methyDesign,sapply(strsplit(as.character(methyDesign$barcode),"-"),function(x) x[6]))
colnames(methyDesign)=c("barcode","platform","patient","subtype","sample","vial","portion","analyte","plate")

methy=do.call(cbind,methySubti)#395741*814
methyRgS=makeGenomicRatioSetFromMatrix(mat=methy,rownames=rownames(methy),
	pData=methyDesign,mergeManifest=T,what="Beta")
methyRgS=dropLociWithSnps(methyRgS,snps="CpG")#polymorphisms may affect DNAm measurements
#394391*814 = filtered 1350 rs from dbSNP
annot=getAnnotation(methyRgS,dropNonMapping=T)#all 394391 probes map to only one locus
methyRgS=methyRgS[annot$chr!="chrY",]#there're 10/13 male LumB & we don't want a fake DM due to chr imbalance 
annot=annot[annot$chr!="chrY",]
methyRgS=methyRgS[annot$chr!="chrX",]
annot=annot[annot$chr!="chrX",]
#384575*814

#Beta values have severe heteroscedasticity for highly methylated or unmethylated CpG sites.
#M-value method provides much better performance in terms of detection rate and true positive rate
#for both highly methylated and unmethylated CpG sites
methySubti=lapply(unique(methyDesign$subtype),function(x) methy[,methyDesign$subtype==x])
methySubti=lapply(methySubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])
names(methySubti)=unique(methyDesign$subtype)
Msubti=lapply(methySubti,beta2m)
sapply(Msubti,dim)
#      Basal   LumA   LumB   Her2 normal
#[1,] 384575 384575 384575 384575 384575
#[2,]    135    331    177     75     96

pdf("methyDistr.pdf")
  plot(density(methySubti$Basal),col="blue",frame.plot=F)
	lines(density(methySubti$LumA),col="purple")
	lines(density(methySubti$LumB),col="green")
	lines(density(methySubti$Her2),col="orange")
	lines(density(methySubti$normal),col="red")
  plot(density(Msubti$Basal),col="blue",frame.plot=F)
  	lines(density(Msubti$LumA),col="purple")
	lines(density(Msubti$LumB),col="green")
	lines(density(Msubti$Her2),col="orange")
	lines(density(Msubti$normal),col="red")
  legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","orange","purple","green","red"),bty="n")
dev.off()
#######correct batches#########################################
library(NOISeq)
noiseqData=readData(data =do.call(cbind,Msubti), factors =methyDesign)
myPCAshared=dat(noiseqData, type = "PCA", norm = T,logtransf=T)
pdf("mPCA.pdf")
 explo.plot(myPCAshared, samples = c(1,2), plottype = "scores", factor = "subtype")
 explo.plot(myPCAshared, samples = c(3,5), plottype = "scores", factor = "plate")
dev.off()

save(methySubti,methyDesign,Msubti,annot,file="ini/imptMethy.RData")
