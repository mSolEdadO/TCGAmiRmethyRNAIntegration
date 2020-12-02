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

#get the data
GDCdownload(xprssn)
expr=GDCprepare(xprssn)
expre=assay(expre)

#keep only samples with metadata
expre=expre[,which(colnames(expre)%in%subtype$samples)]
subtype=subtype[subtype$samples%in%colnames(expre),]
subtype=subtype[order(subtype$patients),]
expre=expre[,order(match(colnames(expre),subtype$samples))]
#there are duplicated measures for 16 patients 
sum(table(subtype[,c(2,4)])>1)
#[1] 16
#keep only transcript id not version numbers
rownames(expre)=sapply(strsplit(rownames(expre),".",fixed=T),
  function(x) x[1])

#table to check for batch effects
designExp=data.frame(cbind(subtype,t(sapply(strsplit(subtype$samples,"-"),function(x) cbind(x)))))
designExp=cbind(designExp,substr(designExp$X4,1,2),
  substr(designExp$X4,3,3))
designExp$X1=designExp$X4=NULL
colnames(designExp)[5:11]=c("TSS","participant","portion_analyte",
  "plate","center","sample","vial")

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
#[1] 19218   867

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
myannot=myannot[myannot$ensembl_gene_id%in%rownames(exprots_hgnc),]

##################CHECK BIASES########################################################
library(NOISeq)
library(edgeR)

#format data for noiseq
noiseqData = readData(data = exprots_hgnc, gc = myannot[,1:2],
 biotype = myannot[,c(1,3)],factor=designExp,
 length=myannot[,c(1,8)])
#1)check expression bias per sample
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#patients with repeated measures
i=colnames(table(subtype[,c(2,4)]))[which(table(subtype[,c(2,4)])>1,
  arr.ind=T)[,2]]
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot",
 samples = which(designExp$patients%in%i))
dev.off()
#2)check for low count genes
myCounts = dat(noiseqData, type = "countsbio", factor = NULL)
png("lowcountsOri.png")
explo.plot(myCounts, plottype = "barplot", 
  samples = which(designExp$patients%in%i))
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
#   761    105 
#explo.plot(mycd,samples=sample(1:ncol(expr),10))

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "center")
png("GCbiasOri.png")
  explo.plot(myGCcontent, samples = NULL)
dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "center")
png("lengthbiasOri.png")
  explo.plot(mylenBias, samples =NULL)
dev.off()
#BUT, since the gene has the same length in all your samples, there is no need to divide by the gene length

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
png("PCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
dev.off()

#################SOLVE BIASES######################################################
#library(cqn)
#library(DESeq2)
#library(pbcmc)
#library("BiocParallel")
#library(sva)
library(EDASeq)

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
#####sigo lo de Diana porque proportion.test, cpm >1 pierde a MIA
countMatrixFiltered = filtered.data(mycountsbio, factor = subtype, norm = FALSE,
 depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#17806 features are to be kept for differential expression analysis with filtering method 1

#2)sigo lo de Diana porque cqn aplana demasiado la varianza
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(designExp,row.names=designExp$samples))
lFull <- withinLaneNormalization(mydataEDA, "length", which = "full")
gcFull <- withinLaneNormalization(lFull, 
  "percentage_gene_gc_content", which = "full")
fullfullTMM <-tmm(normCounts(gcFull), long = 1000, lc = 0, k = 0)

#############################solve batch effect#######################################################
noiseqData = readData(data = fullfullTMM, factors=subtype)
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = T)
ffTMMARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F, norm = "n",  logtransf = T)
myPCA1 = dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = F)
pdf("prepostArsynPCA.pdf")
par(mfrow=c(2,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "subtype")
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "subtype")
dev.off()

#############################final quality check#######################################################
noiseqData = readData(data = exprs(TMMARSyn), factors=subtipos)
mycdnoBatch=dat(noiseqData,type="cd",norm=T)
table(mycdnoBatch@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   496    605 

#final expression dataset
subtipos=subtipos[subtipos$definition!="BRCA.Normal",]
subtipos=as.data.frame(as.matrix(subtipos))
finalData=exprs(TMMARSyn)
finalData=finalData[,colnames(finalData)%in%subtipos$barcode]
save(TMMArsyn,subtipos,file="subtiTMMArsyn.RData")

############################--------------------FIN--------------##########################

#check biotype bias per tissue
#mybiodetection <- dat(noiseqData, k = 0, type = "biodetection", factor = "tissue.definition")
#pdf("ProteinCodingComparison.pdf",height=7,width=7)
#par(oma=c(5,1,1,1))
#explo.plot(mybiodetection,samples=c(1,2),toplot="protein_coding",plottype="comparison")
#[1] "Percentage of protein_coding biotype in each sample:"
#Primary solid Tumor Solid Tissue Normal 
#            35.4141             35.5303 
#[1] "Confidence interval at 95% for the difference of percentages: Primary solid Tumor - Solid Tissue Normal"
#[1] -0.6746  0.4421
#[1] "The percentage of this biotype is NOT significantly different for these two samples (p-value = 0.6868 )."


