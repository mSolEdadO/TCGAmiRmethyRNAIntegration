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
  barcode=subtype$barcode)
GDCdownload(xprssn)
expre=GDCprepare(xprssn,summarizedExperiment=F)#fails if T
temp=as.matrix(expre[,2:ncol(expre)])
rownames(temp)=expre$X1
expre=temp
write.table(expre,"expre.tsv",sep='\t',quote=F)

#subtype to duplicates
i=substr(colnames(expre),1,19)
j=i[duplicated(i)]
designExp=as.matrix(subtype[c(which(!subtype$barcode%in%j),
  as.numeric(sapply(which(subtype$barcode%in%j),rep,2))),])
designExp[designExp[,4]%in%j,1]=colnames(expre)[which(i%in%j)]
designExp=as.data.frame(designExp)
expre=expre[,order(match(colnames(expre),designExp$samples))]
#check all matches

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
#[1] 19218   846

#sum duplicated probes
myannot[myannot$hgnc_id==
  myannot$hgnc_id[duplicated(myannot$hgnc_id)],]
#      ensembl_gene_id percentage_gene_gc_content   gene_biotype start_position
#19876 ENSG00000204510                      50.06 protein_coding       12916610
#55271 ENSG00000279195                      50.06 protein_coding       12916610
#      end_position    hgnc_id hgnc_symbol length
#19876     12920482 HGNC:28415     PRAMEF7   3872
#55271     12920482 HGNC:28415     PRAMEF7   3872
which(rownames(exprots_hgnc)=="ENSG00000279195")
#[1] 19192
exprots_hgnc[rownames(exprots_hgnc)=="ENSG00000204510",]=
  exprots_hgnc[rownames(exprots_hgnc)=="ENSG00000204510",]+
  exprots_hgnc[19192,]
exprots_hgnc=exprots_hgnc[c(1:19191,19193:nrow(exprots_hgnc)),]
dim(exprots_hgnc)
#[1] 19217   846
myannot=myannot[myannot$ensembl_gene_id%in%rownames(exprots_hgnc),]

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
#[1] "Warning: 103 features with 0 counts in all samples are to be removed for this analysis."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   747     98 
#explo.plot(mycd,samples=sample(1:ncol(expr),10))

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
#17359 features are to be kept for differential expression analysis with filtering method 1
myannot=myannot[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered),]

#all names must match
mydataEDA <- newSeqExpressionSet(
  counts=countMatrixFiltered,
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(designExp,row.names=designExp$samples))
lFull <- withinLaneNormalization(mydataEDA, "length", which = "full")#corrects length bias
gcFull <- withinLaneNormalization(lFull, 
  "percentage_gene_gc_content", which = "full")#corrects GC bias & recover length bias
#fullfullTMM <-tmm(normCounts(gcFull), long = 1000, lc = 0, k = 0)
norm.counts <- betweenLaneNormalization(normCounts(gcFull),
 which = "full", offset = FALSE)
#lfull→gcfull→ttm=719 failed
#gcfull→lfull→ttm=731 failed

#############################SOLVE BATCH EFFECT#######################################################
noiseqData = NOISeq::readData(data = norm.counts, factors=designExp)
#has to preceed ARSyN or won't work
mycd=NOISeq::dat(noiseqData,type="cd",norm=TRUE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   200    645 

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
png("Desktop/prepro/lengthbiasFinal.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(mylenBias, samples = x))
dev.off()

#############################RESOLVE DUPLICATES & SAVE##################################################
#get duplicated index
i=designExp$samples[designExp$barcode%in%
  designExp$barcode[duplicated(designExp$barcode)]]
#get sample ids per barcode
i=lapply(designExp$barcode[duplicated(designExp$barcode)],
  function(x) i[substr(i,1,19)%in%x])
#separate duplicates
final=exprs(ffTMMARSyn)
duplis=final[,colnames(final)%in%unlist(i)]
prefi=final[,!colnames(final)%in%unlist(i)]
#average duplicates
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
#identify samples with barcode 
colnames(temp)=designExp$barcode[duplicated(designExp$barcode)]
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
dim(final)
#[1] 17359   842
final=final[,order(match(colnames(final),subtype$barcode))]
write.table(final,"expreNormi.tsv",sep='\t',quote=F)