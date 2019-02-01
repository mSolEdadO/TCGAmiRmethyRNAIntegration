library(TCGAbiolinks)
library(NOISeq)
library(edgeR)
library(DESeq)


#expression matrix: transcripts per sample with the transcript count in each cell
mirna=GDCprepare(mirna)
mirna=mirna[,grep("read_count",colnames(mirna))]
colnames(mirna)=gsub("read_count_","",colnames(mirna))
mirnasN=GDCprepare(mirnasN)
rownames(mirna)=mirnasN$miRNA_ID
mirnasN=mirnasN[,grep("read_count",colnames(mirnasN))]
colnames(mirnasN)=gsub("read_count_","",colnames(mirnasN))
#pimp the table to check for batch effects
mirDesign=cbind(colnames(mirna),substr(colnames(mirna),1,12))
subtipos=subtipos[subtipos$pbcmc2!="normal",]
mirDesign=cbind(mirDesign,sapply(mirDesign[,2],function(x) 
  as.character(subtipos$pbcmc2[subtipos$patient==x])[1]))
mirDesign=rbind(mirDesign,cbind(colnames(mirnasN),substr(colnames(mirnasN),1,12),"normal"))
colnames(mirDesign)=c("barcode","patient","subtype")
mirDesign=as.data.frame(mirDesign)
mirna=cbind(mirna,mirnasN)

noiseqData = readData(data = mirna, factor=mirDesign)
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot",samples=sample(1:ncol(mirna),200))
dev.off()
#check if normalizations is needed
mycd = dat(noiseqData, type = "cd", norm = FALSE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED 
#   880 

#check for low count genes
myCounts = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = sample(1:ncol(mirna),100))
dev.off()
pdf("lowCountThres.pdf")
hist(rowMeans(cpm(mirna,log=T)),ylab="mirna",xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()
#summary(rowSums(mirna>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.0     7.0   103.0   277.3   526.0   881.0 

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
#pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "subtype")
#dev.off()

#filter low counts
#FilteredMatrix = filtered.data(mirna, factor = "subtype",norm = FALSE, method = 3, cpm = 1)
#339 features are to be kept for differential expression analysis with filtering method 3
#con esto me quedo sin miRNAs porque la distribución de lcpm esta muy sesgada a la izq y eso es normal
#it is expected that in miRNA-seq experiments, the 75th percentile of the data will
# be found at only 1 or 2 copies/library [10.1093/bib/bbv019]
sum(rowMedians(as.matrix(mirna))>0)/nrow(mirna)
#[1] 0.2796385 si filtro por mediana de la expresión me quedo con menos de un tercio de la matriz
sum(rowMeans(as.matrix(mirna))>0)/nrow(mirna)
#[1] 0.8442318#lcpmF=cpm(FilteredMatrix,log=T)si filtro por media me quedo casi toda la matriz
#pero puede que los datos sean basura...¿?
FilteredMatrix=mirna[rowMedians(as.matrix(mirna))>0,]
FilteredMatrix1=mirna[rowMeans(as.matrix(mirna))>0,]

#normalize
myTMM=tmm(FilteredMatrix,lc=0)
noiseqData = readData(data = myTMM, factors=mirDesign)
mycdTMM = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   144    736 
myUQUA=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQUA, factors=mirDesign)
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   220    660 
myTMM1=tmm(FilteredMatrix1,lc=0)
noiseqData = readData(data = myTMM1, factors=mirDesign)
mycdTMM1 = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM1@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#880 seguro da bonito porque todo vale 0
myUQUA1=uqua(FilteredMatrix1,lc=0)
noiseqData = readData(data = myUQUA1, factors=mirDesign)
mycdUQUA1 = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA1@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#880 
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix, conditions=mirDesign))
myDESEQ=counts(deseqFactors,normalized=T)
deseqFactors1=estimateSizeFactors(newCountDataSet(FilteredMatrix1,conditions=mirDesign))
myDESEQ1=counts(deseqFactors1,normalized=T)
noiseqData = readData(data = myDESEQ, factors=mirDesign)
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#880 
noiseqData = readData(data = myDESEQ1, factors=mirDesign)
mycdDESEQ1 = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ1@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#880 
noiseqData = readData(data = myTMM1, factors=mirDesign)
TMMARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F, norm = "n",  logtransf = T)

noiseqData = readData(data = exprs(TMMARSyn), factors=mirDesign)
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
explo.plot(mycountsbio, plottype = "boxplot",samples=NULL)#ya se nota la normalización?

myPCA1 = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = TRUE)
pdf("PCA-TMMARSyn.pdf")
par(mfrow = c(2,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,3), plottype = "scores", factor = "subtype")
dev.off()

mirSubti=lapply(unique(mirDesign$subtype),function(x) 
  exprs(TMMARSyn)[,colnames(TMMARSyn)%in%as.character(mirDesign$barcode[mirDesign$subtype==x])])
names(mirSubti)=unique(mirDesign$subtype)
sapply(mirSubti,dim)
#     Basal LumA Her2 LumB normal
#[1,]  1588 1588 1588 1588   1588
#[2,]   142  338  119  178    104
mirSubti=lapply(mirSubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])#keep only one sample per patient
#     Basal LumA Her2 LumB normal
#[1,]  1588 1588 1588 1588   1588
#[2,]   135  331  119  177    104
save(mirSubti,file="ini/mirTMMARSyn.RData")
