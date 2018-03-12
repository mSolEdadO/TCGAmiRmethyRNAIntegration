#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
# &
##    HANDS ON: NGS ANALYSIS - Pre-processing RNA-seq
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

library(biomaRt)  
library(NOISeq)
library(cqn)
library(sva)
library(limma)
library(maSigPro)
library(mdgsa)
library(limma) 
library(edgeR)

#expression matrix: transcripts per sample with the transcript count in each cell
count_files<-list.files("GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Ge
ne_Expression_Quantification/",recursive=T,full.names=T)
count_files <- setNames(count_files, gsub("([:alnum:]+)?.txt", "\\1", basename(c
ount_files)))
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

#annnotate GC content and biotype per transcript
mart=useEnsembl("ensembl",version=77,dataset="hsapiens_gene_ensembl")
#percentage_gene_gc_content for newer version
myannot = getBM(attributes = c("ensembl_gene_id", "percentage_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(count_matrix), mart=mart)
#samples without annotation
#rownames(count_matrix)[!rownames(count_matrix)%in%myannot[,1]]
#keep only annotated transcripts
mygenes = intersect(rownames(count_matrix), myannot$ensembl_gene_id)
count_matrix=count_matrix[rownames(count_matrix)%in%mygenes,]
myannot=myannot[myannot$ensembl_gene_id%in%rownames(count_matrix),]

#calculate gene length=median length for all transcripts
mylength=getBM(attributes=c("ensembl_gene_id","transcript_start","transcript_end"),
	filters="ensembl_gene_id",values=rownames(count_matrix),mart=mart)
transcript_len=abs(mylength[,3]-mylength[,2])
mylength=cbind(mylength,transcript_len)
mylength=sapply(unique(mylength$ensembl_gene_id),function(x) median(mylength[mylength[,1]==x,5]))

#files1=for each sample if it comes from tumor or normal tissue
files1[,2]=unlist(strsplit(files1[,2],'/'))[seq(0,2*nrow(files1),2)]
#to limit tissues to the ones actually present
files1=as.data.frame(as.matrix(files1))

#format data for noiseq
noiseqData = readData(data = count_matrix, gc = myannot[,1:2], biotype = myannot[,c(1,3)],
	factor=files1, length=mylength)

#CHECK BIASES
#1)check expression bias per sample
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = NULL,Colv=NA)
dev.off()
 pdf("noiseqPlot_count_distribution_protein_coding.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = NULL,toplot="protein_coding")
dev.off()

#2)check for low count genes
lcpm=cpm(count_matrix,log=T)
plot(density(lcpm[,colnames(lcpm)%in%files1[files1[,1]=="Primary solid Tumor",2]]),col="cornflowerblue")
lines(density(lcpm[,colnames(lcpm)%in%files1[files1[,1]=="Solid Tissue Normal",2]]),col="green")

#3)check for transcript composition bias
#each sample s is compared to the reference sample r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). Confidence intervals for the M median is
#computed by bootstrapping. If the median of M values for each comparison is not in the CI
# it means that the deviation of the sample is statistically significant. Therefore, a
# normalization procedure should be used 
mycd = dat(noiseqData, type = "cd", norm = FALSE) #slooooow
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "tissue.definition")
pdf("GCbias.pdf")
  explo.plot(myGCcontent, samples = NULL)
 dev.off()
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

###SOLVE BIASES
#1) filter low count genes with the threshold found in 2
countMatrixFiltered = filtered.data(count_matrix, factor = "tissue.definition", 
                       norm = FALSE, method = 3, cpm = 1)
#16143 features are to be kept for differential expression analysis with filtering method 3
lcpmF=cpm(countMatrixFiltered,log=T)
pdf("log-cpm.pdf", width = 5*2)
par(mfrow=c(1,2))
plot(density(lcpm[,colnames(lcpm)%in%files1[files1[,1]=="Primary solid Tumor",2]]),
	col="cornflowerblue",main="raw",xlab="log cpm", bty='n')
lines(density(lcpm[,colnames(lcpm)%in%files1[files1[,1]=="Solid Tissue Normal",2]]),
	col="green")
abline(v=0, lty=3)
legend("topright", legend=unique(files1[,1]), text.col=c("cornflowerblue","green"), bty="n")
plot(density(lcpmF[,colnames(lcpmF)%in%files1[files1[,1]=="Primary solid Tumor",2]]),
	col="cornflowerblue",main="filtered",xlab="log cpm",bty='n')
lines(density(lcpmF[,colnames(lcpmF)%in%files1[files1[,1]=="Solid Tissue Normal",2]]),
	col="green")
abline(v=0, lty=3)
dev.off()

#2)normaliza RNA composition and then GC content coz the 1st needs raw counts
#and gives you scaling factors that can be exported to cqn 
mycdNorm = dat(noiseqData, type = "cd", norm = T)

#3) normalize GC content & gene length
#cqn fits log2(RPM) = s(x) + s(log2(length))
#it needs that all variables have the same ordering
myGC=myannot[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered),1:2]
countMatrixFiltered=countMatrixFiltered[order(match(
	rownames(countMatrixFiltered),names(mylength))),]
mycqn <- cqn(countMatrixFiltered,x = myGC$percentage_gc_content,lengths=1000,
	lengthMethod="fixed",sizeFactors = apply(countMatrixFiltered, 2, sum), verbose = TRUE)
rnaseqCQN = mycqn$y + mycqn$offset   # Note that CQN returns log2-transformed values
rnaseqCQN = rnaseqCQN - min(rnaseqCQN) + 1
noiseqData = readData(data = rnaseqCQN, factors =file, gc = myannot[,c(1,2)],mylength1)
myGCcontentNorm <- dat(noiseqData, k = 0, type = "GCbias", factor = "tissue.definition")
#"Primary solid Tumor":Multiple R-squared:  0.119,  p-value: 0.5099
#"Solid Tissue Normal":Multiple R-squared:  0.1235,  p-value: 0.4754

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
