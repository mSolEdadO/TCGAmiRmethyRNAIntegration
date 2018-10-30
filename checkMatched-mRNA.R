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
library(NOISeq)
library(edgeR)
library(cqn)
library(DESeq2)
library(pbmc)
library("BiocParallel")
library(sva)

#get the data
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary solid Tumor",
  workflow.type = "HTSeq - Counts")
GDCdownload(xprssn)
expr=GDCprepare(xprssn)
xprssN <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = c("Solid Tissue Normal"),
  workflow.type = "HTSeq - Counts")
GDCdownload(xprssN)
normal=GDCprepare(xprssN)
duplicados=colData(expr)$patient[duplicated(colData(expr)$patient)]

#table to check for batch effects
designExp=colData(expr)
designExp=rbind(designExp,colData(normal))
#expression matrix: transcripts per sample with the transcript count in each cell
expr=assay(expr)
expr=cbind(expr,assay(normal))


#check the intersection between normal & tumor samples
venn.diagram(x = list(A=designExp$patient[designExp$definition=="Primary solid Tumor"],
  B=designExp$patient[designExp$definition=="Solid Tissue Normal"]), filename = "Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))

#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "percentage_gene_gc_content", "gene_biotype",
  "start_position","end_position","hgnc_id","hgnc_symbol"),filters = "ensembl_gene_id", 
  values=rownames(expr), mart=mart)
myannot$lenght=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=expr[rownames(expr)%in%myannot$ensembl_gene_id,]
#19236 transcripts

#sum duplicated probes = 2 probes mapping to the same hgnc_id
myannot$hgnc_id[duplicated(myannot$hgnc_id)]
[1] "HGNC:30046"
myannot[myannot$hgnc_id=="HGNC:30046",]
#      ensembl_gene_id percentage_gene_gc_content   gene_biotype start_position
#41743 ENSG00000254093                      43.17 protein_coding       10764963
#44594 ENSG00000258724                      43.86 protein_coding       10725399
#      end_position    hgnc_id hgnc_symbol lenght
#41743     10839884 HGNC:30046       PINX1  74921
#44594     10839847 HGNC:30046       PINX1 114448
exprots_hgnc[rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"],][1,]=
colSums(exprots_hgnc[rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"],])
exprots_hgnc=exprots_hgnc[c(1:(which(rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"])[2]-1),(which(rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"])[2]+1):nrow(exprots_hgnc)),]

##################CHECK BIASES########################################################
#format data for noiseq
noiseqData = readData(data = exprots_hgnc, gc = myannot[,1:2], biotype = myannot[,c(1,3)],
	factor=designExp, length=myannot[,c(1,8)])
#1)check expression bias per sample
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:ncol(expr),200),Colv=NA)
dev.off()
#pdf("noiseqPlot_count_distribution_protein_coding.pdf", width = 15, height = 7)
#explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:889,200),toplot="protein_coding")
#dev.off()
## Biodetection plot
#mybiodetection <- dat(noiseqData, type="biodetection", factor=NULL)
#pdf("mybiodetection.pdf", width = 15, height = 7)
#explo.plot(mybiodetection)
#dev.off()
#saturation plot
mysaturation <- dat(noiseqData, k = 0, ndepth = 7, type = "saturation")
pdf("mysaturation", width = 15, height = 7)
explo.plot(mysaturation, toplot="protein_coding", 
    samples = c(1,12), yleftlim = NULL, yrightlim = NULL)
dev.off()

#2)check for low count genes
myCounts = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = sample(1:ncol(expr),200))
dev.off()
pdf("lowCountThres.pdf")
hist(rowMeans(cpm(exprots_hgnc,log=T)),ylab="genes",xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias
#each sample s is compared to the reference sample r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). Confidence intervals for the M median is
#computed by bootstrapping. If the median of M values for each comparison is not in the CI
# it means that the deviation of the sample is statistically significant. Therefore, a
# normalization procedure should be used 
mycd = dat(noiseqData, type = "cd", norm = FALSE) #slooooow
#[1] "Warning: 739 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#  1081    133 

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "definition")
pdf("GCbias.pdf")
  explo.plot(myGCcontent, samples = NULL)
 dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "definition")
pdf("lengthbias.pdf")
  explo.plot(mylenBias, samples =NULL)
dev.off()
#BUT As the gene has the same length in all your samples, there is no point in dividing by the gene length

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "definition")
dev.off()

#################SOLVE BIASES######################################################
#1) filter low count genes with the threshold found in BIASES(2): The CPM are used because the size 
#of the library per sample affects the power of the experiment. CPM=(counts/fragments
# sequenced)*one million. Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
countMatrixFiltered = filtered.data(exprots_hgnc, factor = "definition", 
norm = FALSE, method = 3, cpm = 1)#<-----------------
#13904 features are to be kept for differential expression analysis with filtering method 3

#2)normaliza RNA composition and then GC content coz the 1st needs raw counts
#and gives you scaling factors that can be exported to the 2nd
factorsTMM=calcNormFactors(countMatrixFiltered)
normTMM=cqn(countMatrixFiltered,lengths=myannot$lenght[
  myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  sizeFactors=factorsTMM,verbose=T)
normalisedTMMatrix=normTMM$y+normTMM$offset
noiseqData = readData(data = normalisedTMMatrix, factors=designExp, gc = myannot[,c(1,2)])
mycdTMM = dat(noiseqData, type = "cd", norm = T)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#cpm=0.1, transcripts=15668 desde aqui me falta MIA (gen de PAM50)
#FAILED PASSED 
#   560    654 
#cpm=1, transcripts=13904 
#FAILED PASSED 
#   430    784 
#cpm=2, transcripts=13164 
#FAILED PASSED 
#   470    744
#cpm=5, transcripts=11908
#FAILED PASSED 
#   476    738
factorsUQUA=calcNormFactors(countMatrixFiltered, method="upperquartile")
normUQUA=cqn(countMatrixFiltered,lengths=myannot$lenght[
  myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  sizeFactors=factorsUQUA,verbose=T)
normalisedUQUAatrix=normUQUA$y+normUQUA$offset
noiseqData = readData(data = normalisedUQUAatrix, factors=designExp, gc = myannot[,c(1,2)])
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   465    749 

rownames(designExp)=designExp$cases
myDESEQ=DESeqDataSetFromMatrix(countData=countMatrixFiltered,colData=designExp,design=~definition)
deseqFactors=estimateSizeFactors(myDESEQ)
normDESEQ=cqn(countMatrixFiltered,lengths=myannot$lenght[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  sizeFactors=sizeFactors(deseqFactors),verbose=T)
myDESEQ=normDESEQ$y+normDESEQ$offset
noiseqData = readData(data = myDESEQ, factors=designExp,gc = myannot[,c(1,2)])
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#cpm=1, transcripts=13904 
#FAILED PASSED 
#   450    764 

#summary(as.numeric(normalizedDESEQmatrix))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  17.33   28.56   30.14   29.84   31.39   38.04 
# summary(as.numeric(normalisedTMMatrix))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  17.21   28.62   30.18   29.89   31.44   38.46 
# summary(as.numeric(normalisedTMMatrix1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  17.23   28.62   30.18   29.89   31.44   38.44 

TMM.TP=normalisedTMMatrix[,colnames(normalisedTMMatrix)%in%designExp$barcode[designExp$shortLetterCode=="TP"]]
TMM.NT=normalisedTMMatrix[,colnames(normalisedTMMatrix)%in%designExp$barcode[designExp$shortLetterCode=="NT"]]

#######################starting classification########################################################
subannot=getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene"),values=rownames(TMM.TP), mart=mart)
colnames(subannot)=c("probe","NCBI.gene.symbol","EntrezGene.ID")#this colnames are needed
subannot=subannot[subannot$probe%in%rownames(TMM.TP),]
subannot=subannot[!duplicated(subannot$probe),]

#49/50 probes are used for clustering
subtypes=molecular.subtyping(sbt.model="pam50",data=t(log2(TMM.TP)),annot=subannot,do.mapping=T)
subtipos=TCGA_MolecularSubtype(colnames(TMM.TP))
subtipos=subtipos$subtypes[,c(4,2)]
temp=cbind(colnames(TMM.TP)[!colnames(TMM.TP)%in%subtipos$barcodes],NA)
colnames(temp)=colnames(subtipos)
subtipos=rbind(subtipos,temp)
subtipos=cbind(subtipos,subtypes$subtype[order(match(names(subtypes$subtype),subtipos$barcode))])
table(subtipos[,2:3],useNA="ifany") 
#                  X
#subtype_PAM50.mRNA Basal Her2 LumB LumA Normal
#     Basal-like       97    0    2    2      0    0
#     HER2-enriched     1   52    5    0      0    0
#     Luminal A         2    6   42  180      2    0
#     Luminal B         1    8  112    4      0    0
#     Normal-like       4    0    2    0      2    0
#     <NA>            100   55  152  251     20    0

#############################solve batch effect#######################################################
temp=cbind(designExp$barcode[!designExp$barcode%in%subtipos[,1]],"normal","normal","normal","normal")
subtipos=rbind(subtipos,temp)

noiseqData = readData(data = normalisedTMMatrix, factors=subtipos, gc = myannot[,c(1,2)])
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = FALSE)
TMMARSyn=ARSyNseq(noiseqData, factor = "genefu", batch = F, norm = "n",  logtransf = FALSE)
myPCA1 = dat(TMM.ARSyn, type = "PCA", norm = T)
pdf("pre-postArsynPCA.pdf")
par(mfrow=c(2,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "genefu")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "genefu")
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "genefu")
explo.plot(myPCA1, samples = c(1,3), plottype = "scores", factor = "genefu")
dev.off()

TMMArsyn.NT=exprs(TP.TMM.ARSyn)[,colnames(exprs(TP.TMM.ARSyn))%in%subtipos[subtipos[,3]=="normal",1]]
TMMArsyn.TP=exprs(TP.TMM.ARSyn)[,colnames(exprs(TP.TMM.ARSyn))%in%subtipos[subtipos[,3]!="normal",1]]

##############final classification based on genefu's alredy classified samples##########################
rownames(subannot)=subannot$probe
obje<-PAM50(exprs=log2(TMMArsyn.TP),annotation=subannot)
obje=filtrate(obje,verbose=T)
obje=classify(obje,verbose=T,std="scale")
obje=permutate(obje,verbose=T,keep=T,nPerm=10000, pCutoff=0.01, where="fdr",corCutoff=0.1,seed=1234567890,
  BPPARAM=MulticoreParam(progressbar=TRUE))
#49/50 probes are used for clustering
temp=permutation(obje)$subtype
subtipos[subtipos[,3]!="normal",5]=temp$Subtype[order(match(rownames(temp),subtipos[,1]))]
#table(subtipos[,c(2,5)],useNA="always")
#               pbcmc2
#tcga            Ambiguous Basal Her2 LumA LumB Normal <NA>
#  Basal-like            0    97    0    2    1      0    1
#  HER2-enriched         0     1   52    0    5      0    0
#  Luminal A             1     2    6  180   37      2    4
#  Luminal B             2     1    8    4   93      0   17
#  Normal-like           0     4    0    0    2      2    0
#  <NA>                  6   100   54  250  131     19   18

--------------noiseqData = readData(data = exprs(TP.TMM.ARSyn), factors=subtipos, gc = myannot[,c(1,2)])
mycdnoBatch=dat(noiseqData,type="cd",norm=T)
table(mycdComBat@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   351    537
write.table(rnaseqCombat[,colnames(rnaseqCombat)%in%rownames(designExp)[designExp$tissue.definition=="Primary solid Tumor"]],
  "GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/mproteTumor.mtrx",sep='\t',quote=F)
write.table(rnaseqCombat[,colnames(rnaseqCombat)%in%rownames(designExp)[designExp$tissue.definition=="Solid Tissue Normal"]],
  "GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/mproteNormal.mtrx",sep='\t',quote=F)




y_DESeq<-DESeqDataSetFromMatrix(countData=rnaseqCombat,colData=designExp, design=~tissue.definition)
cqnOffset <- normTMM$glm.offset
cqnNormFactors <- exp(cqnOffset)
cqnNormFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
normalizationFactors(y_DESeq) <- cqnNormFactors
y_DESeq<-estimateDispersions(y_DESeq, fitType="local")
et_DESeq<-nbinomWaldTest(y_DESeq)
DESeqDEResults<-results(et_DESeq, alpha=0.001)

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


