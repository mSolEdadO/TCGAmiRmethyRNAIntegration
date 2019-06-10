library(TCGAbiolinks)
library(NOISeq)
library(edgeR)
library(DESeq)


#get the data
mirnas <- GDCquery(project = "TCGA-BRCA",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification",
                   barcode=substr(subtipos$barcode,1,12))
GDCdownload(mirnas)
mir=GDCprepare(mirnas)
rownames(mir)=mir$miRNA_ID
mir=mir[,grep("read_count",colnames(mir))]
colnames(mir)=gsub("read_count_","",colnames(mir))

#assign subtype
mirnas=getResults(mirnas)
mirDesign=mirnas[,c(8,27)]
mirDesign=cbind(mirDesign,substr(mirDesign$cases,1,12))
subtipos=subtipos[subtipos$definition!="normal",]
subtipos=cbind(subtipos,substr(subtipos$barcode,1,12))
colnames(mirDesign)[3]="patient"
colnames(subtipos)[3]="patient"
#drop unclassified tumor samples
mirDesign=mirDesign[mirDesign$patient%in%subtipos$patient|mirDesign$tissue.definition=="Solid Tissue Normal",]
mirDesignTumor=mirDesign[mirDesign$tissue.definition!="Solid Tissue Normal",]
mirDesignNormal=mirDesign[mirDesign$tissue.definition=="Solid Tissue Normal",]
temp=sapply(as.character(mirDesignTumor$patient),function(x) as.character(subtipos$definition)[as.character(subtipos$patient)==x])
mirDesignTumor$subtype=temp
mirDesignNormal$subtype="normal"
mirDesignTumor=mirDesignTumor[!duplicated(mirDesignTumor$patient),]
mirDesign=rbind(mirDesignNormal,mirDesignTumor)
mir=mir[,colnames(mir)%in%mirDesign$cases]
mir=mir[,order(match(colnames(mir),mirDesign$cases))]
table(mirDesign$subtype)
#BRCA.Basal  BRCA.Her2  BRCA.LumA  BRCA.LumB     normal 
#       180         81        557        208        104 

#starting quality            
noiseqData = readData(data = mir, factor=mirDesign)
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")#check low counts
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)#check batches
pdf("mirQuality0.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "barplot",samples=1:5)
explo.plot(mycountsbio, plottype = "boxplot",samples=1:5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "subtype")
i=sample(1:ncol(mir),10)
explo.plot(mycd,samples=i)
dev.off()
mycd = dat(noiseqData, type = "cd", norm = FALSE)#check if normalizations is needed
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Warning: 282 features with 0 counts in all samples are to be removed for this analysis."
#FAILED 
#  1129 

#filter low counts
FilteredMatrix = filtered.data(mirna, factor = "subtype",norm = FALSE, method = 1, cpm = 1)
#315 features are to be kept for differential expression analysis with filtering method 3
#con esto me quedo sin miRNAs porque la distribución de lcpm esta muy sesgada a la izq y eso es normal
#it is expected that in miRNA-seq experiments, the 75th percentile of the data will
# be found at only 1 or 2 copies/library [10.1093/bib/bbv019]
temp=lapply(unique(mirDesign$subtype),function(x) 
  mir[,colnames(mir)%in%mirDesign$cases[mirDesign$subtype==x]])
temp1=unique(unlist(lapply(temp,function(x) 
  which(rowSums(x>=5)>=ncol(x)*0.25))))
FilteredMatrix=mir[temp1,]# a minimum of  5  counts  in  at  least  25%  of  the  samples Drago-García2017
#dim(FilteredMatrix1)
#[1] 435 1130

pdf("miRlowCountThres.pdf")
plot(density(rowMeans(cpm(mir[,mirDesign$subtype=="BRCA.Basal"],log=T))),main="no filter",xlab="mean of log CPM",col="blue",ylim=c(0,2.5))
 lines(density(rowMeans(cpm(mir[,mirDesign$subtype=="BRCA.Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(mir[,mirDesign$subtype=="BRCA.LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(mir[,mirDesign$subtype=="BRCA.LumB"],log=T))),col="red")
 lines(density(rowMeans(cpm(mir[,mirDesign$subtype=="normal"],log=T))),col="purple")
legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","green","orange","red","purple"),bty="n")
plot(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="BRCA.Basal"],log=T))),main="aveExpress per condition > 1cpm",xlab="mean of log CPM",col="blue",ylim=c(0,0.15))
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="BRCA.Her2"],log=T))),col="green")
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="BRCA.LumA"],log=T))),col="orange")
 lines(density(rowMeans(cpm(FilteredMatrix[,mirDesign$subtype=="BRCA.LumB"],log=T))),col="red")
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
#    25   1104 
myUQUA=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQUA, factors=mirDesign)
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   401    728 
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix, conditions=mirDesign))
myDESEQ=counts(deseqFactors,normalized=T)
noiseqData = readData(data = myDESEQ, factors=mirDesign)
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    42   1087 
                           
#choose the best normi
noiseqData = readData(data = myTMM, factors=mirDesign)
tmmARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F, norm = "n",  logtransf=T)

#final quality check                         
noiseqData = readData(data = exprs(tmmARSyn), factors=mirDesign)
mycountsbio1 = dat(noiseqData, type = "countsbio", factor = "subtype")
myPCA1 = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = T)
pdf("miRfinalQuali.pdf")
explo.plot(mycdTMM,samples=i)
explo.plot(mycountsbio1, plottype = "boxplot",samples=1:5)#ya se nota la normalización?
explo.plot(mycountsbio1, plottype = "barplot",samples=1:5)#ya se nota la normalización?
par(mfrow = c(1,2))
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,3), plottype = "scores", factor = "subtype")
dev.off()
                           
miRfinal=exprs(tmmARSyn)
colnames(miRfinal)=paste(colnames(miRfinal),mirDesign$subtype,sep=".")
write.table(miRfinal,"normiARSyNmiR.tsv",sep='\t',quote=F)
colnames(mir)=paste(colnames(mir),mirDesign$subtype,sep=".")
write.table(mir,"miR0.tsv",sep='\t',quote=F)

                           

