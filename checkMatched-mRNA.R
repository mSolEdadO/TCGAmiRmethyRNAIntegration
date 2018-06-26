#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
## &&
## Data preparation: 
##      -Quality Control & bias removal
## By Cristobal Fresno - cristobalfresno@gmail.com

library(biomaRt)  
library(NOISeq)
library(cqn)
library(sva)
library(limma)
library(maSigPro)
library(mdgsa)
library(limma) 
library(edgeR)
library(DESeq2)


#expression matrix: transcripts per sample with the transcript count in each cell
count_files<-list.files("GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification",
  recursive=T,full.names=T)
count_files <- setNames(count_files, gsub("([:alnum:]+)?.txt", "\\1", 
  basename(count_files)))
# Initialize outcome variable
count_matrix <- NULL
# Merge all files
for (htseq_file in names(count_files)) {
  file_contents <- read.delim(count_files[htseq_file], "\t",header = FALSE, comment.char = "_")
  colnames(file_contents)[2] <- htseq_file
  # Append the new data to the final matrix
  # Use merge in case there are certain non-shared genes between samples
  count_matrix <- if (is.null(count_matrix)) file_contents else merge(count_matrix, 
                                                                      file_contents,
                                                                      all = TRUE,
																	  by = 1)
}
rownames(count_matrix)=unlist(strsplit(as.character(count_matrix[,1]),'.',fixed=T))[seq(1,2*nrow(count_matrix),2)]
count_matrix[, 1] <- NULL
expresion=getResults(expresion)[,c(3,7,27)]#expresion es el output de GDCquery
expresion=expresion[order(match(expresion$file_name,colnames(count_matrix))),]
colnames(count_matrix)=expresion$cases

#pimp the table to check for batch effects
designExp=expresion[,2:3]
designExp$patients=substr(designExp[,1],1,12)
designExp$sample=substr(designExp[,1],14,15)
designExp$vial=substr(designExp[,1],16,16)
designExp$portion=substr(designExp[,1],18,19)
designExp$analyte=substr(designExp[,1],19,19)
designExp$plate=substr(designExp[,1],22,25)
designExp$center=substr(designExp[,1],27,28)

#check the intersection between normal & tumor samples
venn.diagram(x = list(A=designExp$patients[designExp$tissue.definition=="Primary solid Tumor"],
  B=designExp$patients[designExp$tissue.definition=="Solid Tissue Normal"]), filename = "Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))

#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "percentage_gene_gc_content", "gene_biotype",
  "start_position","end_position","hgnc_id","hgnc_symbol"),filters = "ensembl_gene_id", 
  values=rownames(count_matrix), mart=mart)
myannot$lenght=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]#los duplicados tenian el mismo largo, %gc, inicio y fin, pero distintos hgnc_symbol y/o id 
count_matrix=count_matrix[rownames(count_matrix)%in%myannot$ensembl_gene_id,]

#format data for noiseq
noiseqData = readData(data = count_matrix, gc = myannot[,1:2], biotype = myannot[,c(1,3)],
	factor=designExp, length=myannot[,c(1,8)])

##################CHECK BIASES########################################################
#1)check expression bias per sample
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:889,200),Colv=NA)
dev.off()
 pdf("noiseqPlot_count_distribution_protein_coding.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:889,200),toplot="protein_coding")
dev.off()
## Biodetection plot
mybiodetection <- dat(noiseqData, type="biodetection", factor=NULL)
pdf("mybiodetection.pdf", width = 15, height = 7)
explo.plot(mybiodetection)
dev.off()
#saturation plot
mysaturation <- dat(noiseqData, k = 0, ndepth = 7, type = "saturation")
pdf("mysaturation", width = 15, height = 7)
explo.plot(mysaturation, toplot="protein_coding", 
    samples = c(1,12), yleftlim = NULL, yrightlim = NULL)
dev.off()

#2)check for low count genes
#filter non protein-coding
proteCounts=count_matrix[rownames(count_matrix)%in%
                         myannot$ensembl_gene_id[myannot$gene_biotype=="protein_coding"],]
noiseqData = readData(data = proteCounts, factors=designExp, gc = myannot[,c(1,2)],
  biotype = myannot[,c(1,3)], length=myannot[,c(1,8)])

myCounts = dat(noiseqData, type = "countsbio", factor = "center")
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = 1:100)
dev.off()
pdf("lowCountThres.pdf")
hist(rowMeans(cpm(count_matrix,log=T)),ylab="genes",xlab="mean of log CPM",col="gray")
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
#   739    149 

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "tissue.definition")
pdf("GCbias.pdf")
  explo.plot(myGCcontent, samples = NULL)
 dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "tissue.definition")
pdf("lengthbias.pdf")
  explo.plot(mylenBias, samples =NULL)
dev.off()
#BUT As the gene has the same length in all your samples, there is no point in dividing by the gene length

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "tissue.definition")
dev.off()

#################SOLVE BIASES######################################################
#1) filter low count genes with the threshold found in BIASES(2): The CPM are used because the size 
#of the library per sample affects the power of the experiment. CPM=(counts/fragments
# sequenced)*one million. Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
countMatrixFiltered = filtered.data(proteCounts, factor = "tissue.definition", 
                       norm = FALSE, method = 3, cpm = 1)
#13937 features are to be kept for differential expression analysis with filtering method 3

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
#FAILED PASSED 
#   368    520
#cpm=5, transcripts=12041 
#FAILED PASSED 
#   377    511 

factorsUQUA=calcNormFactors(countMatrixFiltered, method="upperquartile")
normUQUA=cqn(countMatrixFiltered1,lengths=myannot$lenght[
  myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
  sizeFactors=factorsUQUA,verbose=T)
normalisedUQUAatrix=normUQUA$y+normUQUA$offset
noiseqData = readData(data = myUQUA, factors=designExp, gc = myannot[,c(1,2)])
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   765    123
#cpm=5
#FAILED PASSED 
#   774    114 
rownames(designExp)=designExp$cases
#sobreescribiste en cpm=1 con cpm=5, myDESEQ tiene lo de cpm=5
myDESEQ=DESeqDataSetFromMatrix(countData=countMatrixFiltered1,colData=designExp,design=~tissue.definition)
deseqFactors=estimateSizeFactors(myDESEQ)
normDESEQ=cqn(countMatrixFiltered1,lengths=myannot$lenght[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered1)],
  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered1)],
  sizeFactors=sizeFactors(deseqFactors),verbose=T)
myDESEQ=normDESEQ$y+normDESEQ$offset
noiseqData = readData(data = myDESEQ, factors=designExp,gc = myannot[,c(1,2)])
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."
#FAILED PASSED 
#   360    528 
#cpm=5
#FAILED PASSED 
#   373    515 

#summary(as.numeric(normalizedDESEQmatrix))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  18.87   29.14   30.70   30.42   31.94   39.00 
# summary(as.numeric(normalisedTMMatrix))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  18.89   29.17   30.73   30.45   31.97   39.03 
# summary(as.numeric(normalisedTMMatrix1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  18.60   29.89   31.08   31.00   32.18   39.62 
#como las diferencias son mínimas, tomo TMM cpm=1

designCombat = model.matrix(~1,data=designExp)
rnaseqCombat = ComBat(normalisedTMMatrix, batch = designExp$tissue.definition,mod=designCombat,
  par.prior=T) 
noiseqData = readData(data = rnaseqCombat, factors=designExp, gc = myannot[,c(1,2)],length=myannot[,c(1,8)])
myPCAnoBatch = dat(noiseqData, type = "PCA", norm = T)
pdf("noiseqPlot_PCA_after_combat.pdf", width = 5*2, height = 5)
par(mfrow=c(1,2))
explo.plot(myPCAnoBatch, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCAnoBatch, samples = c(1,3), plottype = "scores", factor = "tissue.definition")
dev.off()
mycdComBat=dat(noiseqData,type="cd",norm=T)
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
