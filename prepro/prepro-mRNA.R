#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
## &&
## Data preparation: 
##      -Quality Control & bias removal
## By Cristobal Fresno - cristobalfresno@gmail.com

library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)  
subtype=read.table("subtype.tsv",header=T)

#get the data
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  barcode=subtype$samples)
GDCdownload(xprssn)
expre=GDCprepare(xprssn,summarizedExperiment=F)#fails if T
temp=as.matrix(expre[,2:ncol(expre)])
rownames(temp)=expre$X1
expre=temp
write.table(expre,"RNAseq.tsv",sep='\t',quote=F)

#subtype to duplicates
i=substr(colnames(expre),1,19)
j=i[duplicated(i)]
designExp=subtype[c(which(!subtype$samples%in%j),
  as.numeric(sapply(which(subtype$samples%in%j),rep,2))),]
designExp=designExp[order(match(designExp$samples,substr(colnames(expre),1,19))),]
designExp$barcode=colnames(expre)
#check everything matches

#keep only transcript id not version numbers
rownames(expre)=sapply(strsplit(rownames(expre),".",fixed=T),
  function(x) x[1])

#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
  "percentage_gene_gc_content", "gene_biotype",
  "start_position","end_position","hgnc_id","hgnc_symbol"),
  filters = "ensembl_gene_id", 
  values=rownames(expre),mart=mart)
myannot$length=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&
  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=expre[rownames(expre)%in%myannot$ensembl_gene_id,]
dim(exprots_hgnc)
#[1] 19210   809

#check duplicated probes
#myannot[myannot$hgnc_id==
#  myannot$hgnc_id[duplicated(myannot$hgnc_id)],]
#which(rownames(exprots_hgnc)=="ENSG00000279195")#not anymore
#exprots_hgnc[rownames(exprots_hgnc)=="ENSG00000204510",]=
#  exprots_hgnc[rownames(exprots_hgnc)=="ENSG00000204510",]+
#  exprots_hgnc[19192,]
#exprots_hgnc=exprots_hgnc[c(1:19191,19193:nrow(exprots_hgnc)),]
#dim(exprots_hgnc)
#myannot=myannot[myannot$ensembl_gene_id%in%rownames(exprots_hgnc),]

##################CHECK BIASES########################################################
library(NOISeq)
library(edgeR)

#format data for noiseq
noiseqData = readData(data = exprots_hgnc, gc = myannot[,1:2],
 biotype = myannot[,c(1,3)],factor=designExp,
 length=myannot[,c(1,8)])
#1)check expression bias per subtype
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")
#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()
#2)check for low count genes
png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:5)
dev.off()
png("lowCountThres.png")
hist(rowMeans(cpm(exprots_hgnc,log=T)),ylab="genes",
  xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias
#each sample s is compared to a reference r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). 
#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
mycd = dat(noiseqData, type = "cd", norm = FALSE) #slooooow
#[1] "Warning: 107 features with 0 counts in all samples are to be removed for this analysis."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   695    113 
png("MvaluesOri.png")
 explo.plot(mycd,samples=sample(1:ncol(exprots_hgnc),10))
dev.off()

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias",
 factor = "subtype")
png("GCbiasOri.png",width=1000)
par(mfrow=c(1,5))
 sapply(1:5,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias",
 factor = "subtype")
png("lengthbiasOri.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(mylenBias, samples = x))
dev.off()
#BUT, since the gene has the same length in all your samples, there is no need to divide by the gene length

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = F, logtransf = F)
png("PCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
 factor = "subtype")
dev.off()

#################SOLVE BIASES######################################################
library(EDASeq)

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
#####sigo lo de Diana porque proportion.test, cpm >1 pierde a MIA
countMatrixFiltered = filtered.data(exprots_hgnc, factor = "subtype",
 norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#17077 features are to be kept for differential expression analysis with filtering method 1
myannot=myannot[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered),]

#all names must match
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(designExp,row.names=designExp$barcode))
#order for less bias
gcFull <- withinLaneNormalization(mydataEDA, 
  "percentage_gene_gc_content", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 
fullfullTMM <-NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)
#FAILED PASSED 
#   290    518
noiseqData = NOISeq::readData(data = fullfullTMM,, factors=designExp)
#cd has to preceed ARSyN or won't work
mycd=NOISeq::dat(noiseqData,type="cd",norm=TRUE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   297    511 
#full
#FAILED PASSED 
#   346    462 
#############################SOLVE BATCH EFFECT#######################################################
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = F)
png("preArsyn.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
 factor = "subtype")
dev.off()
ffTMMARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F,
 norm = "n",  logtransf = T)
myPCA = dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = T)
png("postArsyn.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
  factor = "subtype")
dev.off()

#############################FINAL QUALITY CHECK#######################################################
noiseqData = readData(data = exprs(ffTMMARSyn), gc = myannot[,1:2],
 biotype = myannot[,c(1,3)],factor=designExp,
 length=myannot[,c(1,8)])
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype",
  norm=T)
png("CountsFinal.png")
explo.plot(mycountsbio, plottype = "boxplot",samples=1:5)
dev.off()
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", 
  factor = "subtype",norm=T)
png("GCbiasFinal.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", 
  factor = "subtype",norm=T)
png("lengthbiasFinal.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(mylenBias, samples = x))
dev.off()

#############################RESOLVE DUPLICATES & SAVE##################################################
#get duplicates
i=designExp$samples[duplicated(designExp$samples)]
#get sample barcode per sample
i=lapply(i,function(x) designExp$barcode[designExp$samples==x])
#separate duplicates
final=exprs(ffTMMARSyn)
duplis=final[,colnames(final)%in%unlist(i)]
prefi=final[,!colnames(final)%in%unlist(i)]
#average duplicates
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
#identify samples with barcode 
colnames(temp)=designExp$samples[duplicated(designExp$samples)]
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
dim(final)
#[1] 17077   805
final=final[,order(match(colnames(final),subtype$samples))]
write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)

###################################################
                  samples      patient                      barcode
143   TCGA-A7-A13D-01A-13 TCGA-A7-A13D TCGA-A7-A13D-01A-13R-A12P-07
143.1 TCGA-A7-A13D-01A-13 TCGA-A7-A13D TCGA-A7-A13D-01A-13R-A277-07
213   TCGA-A7-A26E-01A-11 TCGA-A7-A26E TCGA-A7-A26E-01A-11R-A277-07
213.1 TCGA-A7-A26E-01A-11 TCGA-A7-A26E TCGA-A7-A26E-01A-11R-A169-07
243   TCGA-A7-A26J-01A-11 TCGA-A7-A26J TCGA-A7-A26J-01A-11R-A169-07
243.1 TCGA-A7-A26J-01A-11 TCGA-A7-A26J TCGA-A7-A26J-01A-11R-A277-07
384   TCGA-A7-A13E-01A-11 TCGA-A7-A13E TCGA-A7-A13E-01A-11R-A12P-07
384.1 TCGA-A7-A13E-01A-11 TCGA-A7-A13E TCGA-A7-A13E-01A-11R-A277-07
########################################################