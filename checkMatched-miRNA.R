library(NOISeq)
library(edgeR)
library(DESeq)


#expression matrix: transcripts per sample with the transcript count in each cell
count_files<-list.files(".",recursive=T,full.names=T)
count_files <- setNames(count_files, gsub("([:alnum:]+)?.txt", "\\1", 
	basename(count_files)))
# Initialize outcome variable
count_matrix <- NULL
# Merge all files
for (htseq_file in names(count_files)) {
  file_contents <- read.delim(count_files[htseq_file], "\t",header = FALSE, 
  	comment.char = "_",colClasses=c("character","integer","NULL","NULL"))
  colnames(file_contents)[2] <- htseq_file
  # Append the new data to the final matrix
  # Use merge in case there are certain non-shared genes between samples
  count_matrix <- if (is.null(count_matrix)) file_contents else merge(count_matrix, 
                                                                      file_contents,
                                                                      all = TRUE,
																	  by = 1)
}
rownames(count_matrix)=count_matrix[,1]
count_matrix[,1]=NULL
count_matrix=count_matrix[rowSums(is.na(count_matrix))<ncol(count_matrix),]
count_matrix=count_matrix[rowSums(count_matrix)>0,]

#pimp the table to check for batch effects
mirDesign=getResults(mirnas)[,c(4,8,29)]#mirnas is the product of TCGAquery
mirDesign$file_name=gsub('.txt','',mirDesign$file_name)
mirDesign$patients=substr(mirDesign$cases,1,12)
mirDesign$sample=substr(mirDesign$cases,14,15)
mirDesign$vial=substr(mirDesign$cases,16,16)
mirDesign$portion=substr(mirDesign$cases,18,19)
mirDesign$analyte=substr(mirDesign$cases,19,19)
mirDesign$plate=substr(mirDesign$cases,22,25)
mirDesign$center=substr(mirDesign$cases,27,28)

noiseqData = readData(data = count_matrix, factor=mirDesign)
#mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#explo.plot(mycountsbio, plottype = "boxplot",samples=NULL)

#check if normalizations is needed
mycd = dat(noiseqData, type = "cd", norm = FALSE)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."

#check for low count genes
lcpm=cpm(count_matrix, log=T)
plot(density(lcpm[,colnames(lcpm)%in%mirDesign[mirDesign[,2]=="Primary solid Tumor",3]]),col="cornflowerblue")
lines(density(lcpm[,colnames(lcpm)%in%mirDesign[mirDesign[,2]=="Solid Tissue Normal",3]]),col="green")
#sum(count_matrix==0)/length(as.matrix(count_matrix))#hay muchos mirnas sin cuentas

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "tissue.definition")
dev.off()

#filter low counts
#FilteredMatrix = filtered.data(count_matrix, factor = "tissue.definition", 
#                       norm = FALSE, method = 3, cpm = 10)
#con esto temo quedarme sin miRNAs porque la distribución de lcpm esta muy sesgada a la izq 
#en su lugar hago como Diana Drago: filter targets with a minimum of  5  counts  in  at  
#least  25%  of  the  samples; pero bajando el min a 1

temp=count_matrix[rowSums(count_matrix>0)>(886*0.25),]#tambien jala con %50 y %75

#lcpmF=cpm(FilteredMatrix,log=T)
#pdf("lowCounts.pdf", width = 5*2)
#par(mfrow=c(1,2))
#plot(density(lcpm[,colnames(lcpm)%in%mirDesign[mirDesign[,2]=="Primary solid Tumor",3]]),
#	col="cornflowerblue",main="raw",xlab="log cpm", bty='n',ylim=c())
#lines(density(lcpm[,colnames(lcpm)%in%mirDesign[mirDesign[,2]=="Solid Tissue Normal",3]]),
#	col="green")
#abline(v=-2.6, lty=3)
#legend("topright", legend=unique(mirDesign[,2]), text.col=c("cornflowerblue","green"), bty="n")
#plot(density(lcpmF[,colnames(lcpmF)%in%mirDesign[mirDesign[,2]=="Primary solid Tumor",3]]),
#	col="cornflowerblue",main="filtered",xlab="log cpm",bty='n')
#lines(density(lcpmF[,colnames(lcpmF)%in%mirDesign[mirDesign[,2]=="Solid Tissue Normal",3]]),
#	col="green")
#abline(v=-2.6, lty=3)
#dev.off()

#normalize
myTMM=tmm(temp,lc=0)
noiseqData = readData(data = myTMM, factors=mirDesign)
mycdTMM = dat(noiseqData, type = "cd", norm = T)
#[1] "Diagnostic test: PASSED."

#deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix,
#	conditions=mirDesign[mirDesign[,3]%in%colnames(FilteredMatrix),c(3,2)]))
#myDESEQ=counts(deseqFactors,normalized=T)
#noiseqData = readData(data = myDESEQ, factors=mirDesign[mirDesign[,3]%in%colnames(FilteredMatrix),c(3,2)])
#mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#explo.plot(mycountsbio, plottype = "boxplot",samples=NULL)#ya se nota la normalización
#mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
#table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#myUQUA=uqua(FilteredMatrix,lc=0)
#noiseqData = readData(data = myUQUA, factors=mirDesign[mirDesign[,3]%in%colnames(FilteredMatrix),c(3,2)])
#mycdUQUA = dat(noiseqData, type = "cd", norm = T)
#temp=dat2save(mycdUQUA)
#temp1=cbind(temp$DiagnosticTest[match(mirDesign[,3],rownames(temp$DiagnosticTest)),3],mirDesign[,2])
#table(as.data.frame(temp1))
#      tumor normal  
#  FAILED 7  14
#  PASSED 91 78


myPCA = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = TRUE)
pdf("PCAnormalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "tissue.definition")
dev.off()
#hace falta corregir batches?

mirDesign=mirDesign[order(match(mirDesign$file_name,colnames(myTMM))),]
colnames(myTMM)=mirDesign$cases
mirTumor=myTMM[,colnames(myTMM)%in%mirDesign$cases[mirDesign$tissue.definition=="Primary solid Tumor"]]
#como el diagnostico dio "PASSED" para todas las muestras, tomo cualquiera para los pacientes repetidos
mirTumor=mirTumor[,colnames(mirTumor)%in%mirDesign$cases[mirDesign$cases%in%colnames(mirTumor)][!duplicated(mirDesign$patients[mirDesign$cases%in%colnames(mirTumor)])]]
mirNormal=myTMM[,colnames(myTMM)%in%mirDesign$cases[mirDesign$tissue.definition=="Solid Tissue Normal"]]
write.table(mirNormal,"mirNormal.mtrx",sep='\t',quote=F)
write.table(mirTumor,"mirTumor.mtrx",sep='\t',quote=F)
