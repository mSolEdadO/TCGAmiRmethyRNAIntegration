library(TCGAbiolinks)
library(SummarizedExperiment)
library(BEclear)
library(NOISeq)
library(sva)


#read matrix
methy=GDCprepare(mthyltn, save=T, save.filename="methy.rda", directory="GDCdata")
#filter probes with NA
load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy.rda")
meth=unlist(assays(data))
meth=meth[rowSums(is.na(meth))==0,]

#pimp the table to check for batch effects
methDesign=list()
methDesign$name=colnames(data)
methDesign$patients=substr(methDesign$name,1,12)
methDesign$sample=substr(methDesign$name,14,15)
methDesign$vial=substr(methDesign$name,16,16)
methDesign$portion=substr(methDesign$name,18,19)
methDesign$analyte=substr(methDesign$name,19,19)
methDesign$plate=substr(methDesign$name,22,25)
methDesign$center=substr(methDesign$name,27,28)
methDesign=as.data.frame(sapply(methDesign,function(x) x))

#check the intersection between normal & tumor samples
venn.diagram(x = list(A=methDesign$patients[methDesign$sample=="01"],
	B=methDesign$patients[methDesign$sample=="11"]), filename = "Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))

#check normalization
normal=as.numeric(meth[,colData(methy)$definition=="Solid Tissue Normal"])
tumor=as.numeric(meth[,colData(methy)$definition=="Primary solid Tumor"])
png("normal?.png")
qqplot(normal,tumor)
dev.off()

#check batches
noiseqData = readData(data = meth, factors = methDesign)
pcaMth = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = T)
pdf("batchMethy.pdf")
explo.plot(pcaMth, samples = c(1,2), plottype = "scores", factor = "sample")
 explo.plot(pcaMth, samples = c(1,3), plottype = "scores", factor = "sample")
dev.off()
#tumor samples are different from normal samples

#
M=log2(meth/(1-meth))
batchModel=model.matrix(~1,data=methDesign)
mthCombat = ComBat(M, batch =methDesign$sample,mod =batchModel, par.prior = F)
noiseqData = readData(data = mthCombat, factors = methDesign)
pcaMth = dat(noiseqData, type = "PCA", norm =T, logtransf = T)
pdf("batchMethyN.pdf")
explo.plot(pcaMth, samples = c(1,2), plottype = "scores", factor = "sample")
 explo.plot(pcaMth, samples = c(1,3), plottype = "scores", factor = "sample")
dev.off()

#divide data in tumor and normal matrixes
temp=methDesign$name[methDesign$sample=="01"]
tumor=mthCombat[,colnames(mthCombat)%in%temp[!duplicated(methDesign$patients[methDesign$sample=="01"])]]
temp=methDesign$name[methDesign$sample=="11"]
write.table(tumor,"methyTumor.mtrx",sep='\t',quote=F)
temp=methDesign$name[methDesign$sample=="11"]
normal=mthCombat[,colnames(mthCombat)%in%temp[!duplicated(methDesign$patients[methDesign$sample=="11"])]]
 write.table(normal,"methyNormal.mtrx",sep='\t',quote=F)

