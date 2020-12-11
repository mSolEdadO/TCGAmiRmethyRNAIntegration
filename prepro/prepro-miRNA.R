library(TCGAbiolinks)

subtype=read.table("Desktop/prepro/subtype.tsv",header=T)
#get the data
mirnas <- GDCquery(project = "TCGA-BRCA",
	data.category = "Transcriptome Profiling",
	data.type = "miRNA Expression Quantification",
	barcode=subtype$barcode)
GDCdownload(mirnas)
mir=GDCprepare(mirnas)
rownames(mir)=mir$miRNA_ID
mir=mir[,grep("read_count",colnames(mir))]
colnames(mir)=gsub("read_count_","",colnames(mir))
dim(mir)
#[1] 1881  846
write.table(mir,"Desktop/prepro/miR.tsv",sep='\t',quote=F)

#subtype to duplicates
i=substr(colnames(mir),1,19)
j=i[duplicated(i)]
designExp=as.matrix(subtype[c(which(!subtype$barcode%in%j),
  as.numeric(sapply(which(subtype$barcode%in%j),rep,2))),])
designExp=as.data.frame(designExp)
designExp=designExp[order(match(designExp$barcode,
	substr(colnames(mir),1,19))),]

##################CHECK BIASES########################################################
library(NOISeq)

noiseqData = readData(data = mir, factor=designExp)
mycountsbio = dat(noiseqData, type = "countsbio",factor = "subtype")#check low counts
#distributions
png("miROri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()
#counts
png("miRcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:5)
dev.off()
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, 
logtransf = FALSE)#check batches
png("miRPCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
	factor = "subtype")
dev.off()
mycd = dat(noiseqData, type = "cd", norm = FALSE)#check if normalizations is needed
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Warning: 295 features with 0 counts in all samples are to be removed for this analysis."
#FAILED 
#   845 
png("miRcdOri.png")
 explo.plot(mycd,samples=1:10)
dev.off()

#################SOLVE BIASES######################################################
library(DESeq)
library(EDAseq)

#filter low counts
FilteredMatrix = filtered.data(mir, factor = "subtype",
	norm = FALSE, method = 1, cpm = 0)
#669 features are to be kept for differential expression analysis with filtering method 1
#it is expected that in miRNA-seq experiments, the 75th percentile 
#of the data will be found at only 1 or 2 copies/library [10.1093/bib/bbv019]
#Drago-García2017 used a minimum of  5  counts  in  at  least  25%
# of  the  samples 
#temp=lapply(unique(designExp$subtype),function(x) 
#  mir[,colnames(mir)%in%designExp$samples[designExp$subtype==x]])
#temp1=unique(unlist(lapply(temp,function(x) 
#  which(rowSums(x>=5)>=ncol(x)*0.25))))
#length(temp1)
#[1] 446

#normalize
#UQ, TMM, DESeq, quantile and cyclic loess normalized data have the
#most similar fold-change estimates
#10.1093/bib/bbv019 
myTMM=tmm(FilteredMatrix,lc=0)
noiseqData = readData(data = myTMM, factors=designExp)
mycdTMM = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   810     35 
myUQUA=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQUA, factors=designExp)
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   828     17 
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix,
 conditions=designExp))
myDESEQ=counts(deseqFactors,normalized=T)
noiseqData = readData(data = myDESEQ, factors=designExp)
mycdDESEQ = NOISeq::dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   733    112 
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(FilteredMatrix),
  phenoData=data.frame(designExp,row.names=designExp$samples))
norm.counts <- betweenLaneNormalization(mydataEDA,
 which = "full", offset = FALSE)
noiseqData = readData(data = assayData(norm.counts)$normalizedCounts,
 factors=designExp)
mycdFull = dat(noiseqData, type = "cd", norm = T)
table(mycdFull@dat$DiagnosticTest[,  "Diagnostic Test"])
#PASSED 
#   845
png("miRcd_final.png")
 explo.plot(mycdFull,samples=1:10)
dev.off()
#############################SOLVE BATCH EFFECT#######################################################

nobatch=ARSyNseq(noiseqData, factor = "subtype", batch = F,
 norm = "n",  logtransf=F)

#############################FINAL QUALITY CHECK#######################################################
noiseqData = readData(data = exprs(nobatch), factors=designExp)
mycountsbio = dat(noiseqData, type = "countsbio",factor = "subtype")#check low counts
png("miRFinal.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()
png("miRcountsFinal.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:5)
dev.off()
myPCA = dat(noiseqData, type = "PCA", norm = T,logtransf=F)#check batches
png("miRPCA_Final.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
	factor = "subtype")
dev.off()
#############################RESOLVE DUPLICATES & SAVE##################################################
miRfinal=exprs(nobatch)
#get duplicated index
i=designExp$samples[designExp$barcode%in%
  designExp$barcode[duplicated(designExp$barcode)]]
#get sample ids per barcode
i=lapply(designExp$barcode[duplicated(designExp$barcode)],
  function(x) i[substr(i,1,19)%in%x])
#separate duplicates
duplis=miRfinal[,colnames(miRfinal)%in%unlist(i)]
prefi=miRfinal[,!colnames(miRfinal)%in%unlist(i)]
#average duplicates
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
#identify samples with barcode 
colnames(temp)=designExp$barcode[duplicated(designExp$barcode)]
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
dim(final)
#[1] 669 842
final=final[,order(match(colnames(final),subtype$barcode))]
write.table(final,"miRNormi.tsv",sep='\t',quote=F)                           
                           

