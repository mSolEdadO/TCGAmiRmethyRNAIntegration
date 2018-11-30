library(TCGAbiolinks)
library(SummarizedExperiment)
library(NOISeq)
library(Venn.Diagram)
library(BEclear)
library(data.table)

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
methySubti=sapply(c("Basal","LumA","LumB","Her2"),function(x) 
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
methySubti=lapply(methySubti,function(x) impute.knn(x,k=15)$data)
#Buuren and Groothuis suggested that 15 variables are generally sufficient for imputation

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
shared=do.call(cbind,sapply(methySubti,function(x) x[rownames(x)%in%probes,]))
#mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#pdf("methy_ori.pdf", width = 15, height = 7)
#explo.plot(mycountsbio, plottype = "boxplot", samples = sample(which(methyDesign$shortLetterCode!='NT'),200),Colv=NA)
#explo.plot(mycountsbio, plottype = "boxplot", samples = which(methyDesign$shortLetterCode=='NT'),Colv=NA)
#dev.off()
samples=data.table(methyDesign[,c(1,4)])
bedata=data.table(feature=rownames(shared),shared)
bedata=melt(data=bedata,id.vars="feature",variable.name="sample",value.name="beta.value")
setkey(bedata,"feature","sample")
pvals=calcPvalues(bedata,samples,adjust=T,method=T)# ks test for every gene in every batch
#every gene with a p-value below 0.01 is assumed as batch affected
mdifs <- calcMedians(data, samples)#median difference (mdif) values for every gene and  every  batch
#consider  all  genes  with  mdif  values beyond or equal to 0.05 as batch affected since values beyond
# this threshold would  already  make  a  biological  difference  in  the  methylation  level
summary <- calcSummary(medians=mdifs, pvalues=pvals)
score <- calcScore(ex.data, ex.samples, summary, dir=getwd())#for every batch if it is batch effected
#boxplot para algunas muestras
cleared.data <- clearBEgenes(ex.data, ex.samples, summary)#set affected beta values to NA
corrected.data <- BEclear(cleared.data, rowBlockSize=60, colBlockSize=60,
	epochs=50, outputFormat="RData", dir=getwd())#matrix completion using Latent Factor Models based on
#matrix factorization
#replaceWrongValues() #during  the  prediction, values  beyond  the boundaries of beta values can be returned








myPCA2 = dat(sharedARSyn, type = "PCA", norm = T)
pdf("batchMethy450.pdf")
par(mfrow=c(2,2))
explo.plot(pcaMth, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(pcaMth, samples = c(1,3), plottype = "scores", factor = "subtype")
explo.plot(myPCA2, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA2, samples = c(1,3), plottype = "scores", factor = "subtype")
dev.off()

#normal=methySubti$normal
#Her2=methySubti$Her2

#drop probe annotation not present anymore
annotMethy450=annotMethy[annotMethy[,1]%in%rownames(shared),]
#extract subtype matrix
methySubti=lapply(unique(methyDesign$subtype),function(x) 
	sharedARSyn[,colnames(sharedARSyn)%in%methyDesign$barcode[methyDesign$subtype==x]])
#methySubti$Her2=NULL
methySubti=lapply(methySubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])#keep only one sample per patient
#sapply(methySubti,dim)
#       LumA  Basal   LumB normal
#[1,] 395808 395808 395808 395808
#[2,]    331    135    177     96
save(annotMethy450,methySubti,file="ini/methyImputNoBatch450.RData")

#due to batch correction, data is alredy outside of the [0,1] beta range, and so the next step is unnecesary 
#Beta values have severe heteroscedasticity for highly methylated or unmethylated CpG sites.
#M-value method provides much better performance in terms of detection rate and true positive rate
#for both highly methylated and unmethylated CpG sites
#M=log2((methySubti$Her2+min(methySubti$Her2))/(1-(methySubti$Her2+min(methySubti$Her2))))

#batches on 27K platform data
#Her2=Her2[rownames(Her2)%in%rownames(normal),]
#normal=normal[rownames(normal)%in%rownames(Her2),]
#normal=normal[order(match(rownames(normal),rownames(Her2))),]
#shared=cbind(Her2,normal)#Her2 needs its own reference coz the dataset is different
#temp=as.data.frame(as.matrix(methyDesign[methyDesign$barcode%in%colnames(shared),]))
#temp=cbind(temp,paste(temp$subtype,temp$platform,sep=':'))
#colnames(temp)[5]="subtype:platform"

#noiseqData=readData(data=shared,factor=temp)
#pcaHer2 = dat(noiseqData, type = "PCA", norm = T, logtransf = T)
#designCombat=model.matrix(~temp$subtype)#so I don't lose the subtype signal while correcting for the platform batch
#nobatch=ComBat(shared,batch=temp$platform,mod=designCombat,par.prior=T)
#noiseqData=readData(data=nobatch,factor=temp)
#pcaNobatch = dat(noiseqData, type = "PCA", norm = T, logtransf = T)
#pdf("batchMethyher2.pdf")
#par(mfrow=c(2,2))
#explo.plot(pcaHer2, samples = c(1,2), plottype = "scores", factor = "subtype:platform")
#explo.plot(pcaHer2, samples = c(1,3), plottype = "scores", factor = "subtype:platform")
#explo.plot(pcaNobatch, samples = c(1,2), plottype = "scores", factor = "subtype:platform")
#explo.plot(pcaNobatch, samples = c(1,3), plottype = "scores", factor = "subtype:platform")
#dev.off()

#annotMethy27=annotMethy[annotMethy$Composite.Element.REF%in%rownames(nobatch),]
#methy27=sapply(unique(temp$subtype),function(x) nobatch[,colnames(nobatch)%in%temp$barcode[temp$subtype==x]])
#save(annotMethy27,methy27,methyDesign,file="ini/methyImputNoBatch27.RData")
#sapply(methy27,function(x) sum(duplicated(substr(colnames(x),1,12))))
#  Her2 normal no duplicated patient samples
#     0      0 

