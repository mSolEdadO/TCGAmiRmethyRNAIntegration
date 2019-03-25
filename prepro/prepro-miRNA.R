library(TCGAbiolinks)
library(NOISeq)
library(edgeR)
library(DESeq)


#expression matrix: transcripts per sample with the transcript count in each cell
mirna=GDCprepare(mirna)
rownames(mirna)=mirna[,1]
mirna=mirna[,grep("read_count",colnames(mirna))]
colnames(mirna)=gsub("read_count_","",colnames(mirna))
mirnasN=GDCprepare(mirnasN)
mirnasN=mirnasN[,grep("read_count",colnames(mirnasN))]
colnames(mirnasN)=gsub("read_count_","",colnames(mirnasN))
#pimp the table to check for batch effects
mirDesign=cbind(colnames(mirna),substr(colnames(mirna),1,12))
subtipos=subtipos[subtipos$pbcmc2!="normal",]
mirDesign=cbind(mirDesign,sapply(mirDesign[,2],function(x) 
  as.character(subtipos$pbcmc2[subtipos$patient==x])[1]))
mirDesign=rbind(mirDesign,cbind(colnames(mirnasN),substr(colnames(mirnasN),1,12),"normal"))
colnames(mirDesign)=c("sample","patient","subtype")
mirDesign=as.data.frame(mirDesign)
mirna=cbind(mirna,mirnasN)

noiseqData = readData(data = mirna, factor=mirDesign)
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
#explo.plot(mycountsbio, plottype = "boxplot",samples=sample(1:ncol(mirna),200))
#dev.off()
#check if normalizations is needed
mycd = dat(noiseqData, type = "cd", norm = FALSE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Warning: 298 features with 0 counts in all samples are to be removed for this analysis."
#FAILED 
#   836 

#check for low count genes
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = sample(1:ncol(mirna),100))
dev.off()
#summary(rowSums(mirna>0))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.0     6.0    99.0   264.4   500.0   837.0 

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
#par(mfrow = c(1,2))
#explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
#explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "subtype")
#dev.off()

#filter low counts
FilteredMatrix = filtered.data(mirna, factor = "subtype",norm = FALSE, method = 1, cpm = 1)
#323 features are to be kept for differential expression analysis with filtering method 3
#con esto me quedo sin miRNAs porque la distribución de lcpm esta muy sesgada a la izq y eso es normal
#it is expected that in miRNA-seq experiments, the 75th percentile of the data will
# be found at only 1 or 2 copies/library [10.1093/bib/bbv019]
sum(rowMedians(as.matrix(mirna))>0)/nrow(mirna)
#[1] 0.281765 si filtro por mediana de la expresión me quedo con menos de un tercio de la matriz
sum(rowMeans(as.matrix(mirna))>0)/nrow(mirna)
#[1] 0.8415736#lcpmF=cpm(FilteredMatrix,log=T)si filtro por media me quedo casi toda la matriz
#pero puede que los datos sean basura...¿?
temp=lapply(unique(mirDesign$subtype),function(x) mirna[,colnames(mirna)%in%mirDesign$sample[mirDesign$subtype==x]])
temp1=unique(unlist(lapply(temp,function(x) which(rowSums(x>=5)>=ncol(x)*0.25))))
FilteredMatrix1=mirna[temp1,]# a minimum of  5  counts  in  at  least  25%  of  the  samples Drago-García2017
#dim(FilteredMatrix1)
#[1] 446 837
temp1=unique(unlist(lapply(temp,function(x) which(rowSums(x>=1)>=ncol(x)*0.25))))
FilteredMatrix2=mirna[temp1,]# a minimum of  1  counts  in  at  least  25%  of  the  samples Drago-García2017

pdf("lowCountThres.pdf")
plot(density(rowMeans(cpm(mirna[,mirDesign$subtype=="Basal"],log=T))),main="no filter",xlab="mean of log CPM",col="blue",ylim=c(0,1))
 lines(density(rowMeans(cpm(mirna[,mirDesign$subtype=="Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(mirna[,mirDesign$subtype=="LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(mirna[,mirDesign$subtype=="LumB"],log=T))),col="red")
 lines(density(rowMeans(cpm(mirna[,mirDesign$subtype=="normal"],log=T))),col="purple")
legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","green","orange","red","purple"),bty="n")
plot(density(rowMeans(cpm(FilteredMatrix1[,mirDesign$subtype=="Basal"],log=T))),main="25% per condition > 5",xlab="mean of log CPM",col="blue",ylim=c(0,0.15))
 lines(density(rowMeans(cpm(FilteredMatrix1[,mirDesign$subtype=="Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(FilteredMatrix1[,mirDesign$subtype=="LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(FilteredMatrix1[,mirDesign$subtype=="LumB"],log=T))),col="red")
 lines(density(rowMeans(cpm(FilteredMatrix1[,mirDesign$subtype=="normal"],log=T))),col="purple")
legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","green","orange","red","purple"),bty="n")
plot(density(rowMeans(cpm(FilteredMatrix2[,mirDesign$subtype=="Basal"],log=T))),main="25% per condition > 1",xlab="mean of log CPM",col="blue",ylim=c(0,0.2))
 lines(density(rowMeans(cpm(FilteredMatrix2[,mirDesign$subtype=="Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(FilteredMatrix2[,mirDesign$subtype=="LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(FilteredMatrix2[,mirDesign$subtype=="LumB"],log=T))),col="red")
 lines(density(rowMeans(cpm(FilteredMatrix2[,mirDesign$subtype=="normal"],log=T))),col="purple")
legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","green","orange","red","purple"),bty="n")
plot(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="Basal"],log=T))),main="aveExpress per condition > 1cpm",xlab="mean of log CPM",col="blue",ylim=c(0,0.2))
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="LumB"],log=T))),col="red")
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="normal"],log=T))),col="purple")
legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","green","orange","red","purple"),bty="n")
dev.off()

#normalize
#Adjusting the data by cpm or total count scaling introduces more variability, whereas 
#UQ and TMM decreased the variance across all miRNAs
#Data normalized by UQ or TMM clearly dominated the Precision-recall curves
#UQ, TMM, DESeq, quantile and cyclic loess normalized data have the most similar fold-change estimates
#10.1093/bib/bbv019 
myTMM=tmm(FilteredMatrix,lc=0)
noiseqData = readData(data = myTMM, factors=mirDesign)
mycdTMM = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    15    821 
myUQUA=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQUA, factors=mirDesign)
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   249    587
myTMM1=tmm(FilteredMatrix1,lc=0)
noiseqData = readData(data = myTMM1, factors=mirDesign)
mycdTMM1 = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM1@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    16    820 vs 821, me quedo con ésta
myUQUA1=uqua(FilteredMatrix1,lc=0)
noiseqData = readData(data = myUQUA1, factors=mirDesign)
mycdUQUA1 = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA1@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   168    668 
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix, conditions=mirDesign))
myDESEQ=counts(deseqFactors,normalized=T)
deseqFactors1=estimateSizeFactors(newCountDataSet(FilteredMatrix1,conditions=mirDesign))
myDESEQ1=counts(deseqFactors1,normalized=T)
noiseqData = readData(data = myDESEQ, factors=mirDesign)
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    36    800 
noiseqData = readData(data = myDESEQ1, factors=mirDesign)
mycdDESEQ1 = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ1@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    85    751 

pdf("miRsampleBias.pdf")
temp=sample(1:ncol(mirna),10)
explo.plot(mycd,samples=temp)
explo.plot(mycdTMM1,samples=temp)
dev.off()
noiseqData = readData(data = myTMM1, factors=mirDesign)
tmmARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F, norm = "n",  logtransf =F)

noiseqData = readData(data = exprs(tmmARSyn), factors=mirDesign)
mycountsbio1 = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot",samples=sample(1:ncol(mirna),200))
explo.plot(mycountsbio1, plottype = "boxplot",samples=sample(1:837,200))#ya se nota la normalización?
dev.off()

myPCA1 = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = TRUE)
pdf("PCA-TMMARSyn.pdf")
par(mfrow = c(2,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,3), plottype = "scores", factor = "subtype")
dev.off()

mirSubti=lapply(unique(mirDesign$subtype),function(x) 
  exprs(tmmARSyn)[,colnames(tmmARSyn)%in%as.character(mirDesign$sample[mirDesign$subtype==x])])
names(mirSubti)=unique(mirDesign$subtype)
sapply(mirSubti,dim)
#     LumA Her2 Basal LumB normal
#[1,]  446  446   446  446    446
#[2,]  338   75   142  178    104
mirSubti=lapply(mirSubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])#keep only one sample per patient
#     Basal LumA Her2 LumB normal
#[1,]   446  446 446  446    446
#[2,]   135  331  75  177    104
save(mirSubti,file="ini/mirTMMARSyn.RData")


