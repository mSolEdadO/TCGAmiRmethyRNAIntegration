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

#which is tumor and which is normal
mydesign=read.table("/home/msoledad/Documents/tordoc/data/archivos.tsv",sep='\t',header=T)
#pimp mydesign
mydesign=mydesign[mydesign$experimental_strategy=="miRNA-Seq",]
mydesign$file=unlist(strsplit(as.character(mydesign$file),"/"))[seq(7,7*nrow(mydesign),7)]
mydesign$file=gsub(".txt","",mydesign$file)
mydesign$cases=substr(mydesign$cases,1,12)
mydesign$experimental_strategy=NULL

noiseqData = readData(data = count_matrix, factor=mydesign[,c(3,2)])
#mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#explo.plot(mycountsbio, plottype = "boxplot",samples=NULL)

#check if normalizations is needed
mycd = dat(noiseqData, type = "cd", norm = FALSE)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."

#check for low count genes
lcpm=cpm(count_matrix)
plot(density(lcpm[,colnames(lcpm)%in%mydesign[mydesign[,2]=="Primary solid Tumor",3]]),col="cornflowerblue")
lines(density(lcpm[,colnames(lcpm)%in%mydesign[mydesign[,2]=="Solid Tissue Normal",3]]),col="green")
#sum(count_matrix==0)/length(as.matrix(count_matrix))#hay muchos mirnas sin cuentas

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "tissue.definition")
dev.off()

#filter low counts
FilteredMatrix = filtered.data(count_matrix, factor = "tissue.definition", 
                       norm = FALSE, method = 3, cpm = 10)
#205 features are to be kept for differential expression analysis with filtering method 3
lcpmF=cpm(FilteredMatrix,log=T)
pdf("lowCounts.pdf", width = 5*2)
par(mfrow=c(1,2))
plot(density(lcpm[,colnames(lcpm)%in%mydesign[mydesign[,2]=="Primary solid Tumor",3]]),
	col="cornflowerblue",main="raw",xlab="log cpm", bty='n',ylim=c())
lines(density(lcpm[,colnames(lcpm)%in%mydesign[mydesign[,2]=="Solid Tissue Normal",3]]),
	col="green")
abline(v=-2.6, lty=3)
legend("topright", legend=unique(mydesign[,2]), text.col=c("cornflowerblue","green"), bty="n")
plot(density(lcpmF[,colnames(lcpmF)%in%mydesign[mydesign[,2]=="Primary solid Tumor",3]]),
	col="cornflowerblue",main="filtered",xlab="log cpm",bty='n')
lines(density(lcpmF[,colnames(lcpmF)%in%mydesign[mydesign[,2]=="Solid Tissue Normal",3]]),
	col="green")
abline(v=-2.6, lty=3)
dev.off()

#normalize
myTMM=tmm(FilteredMatrix,lc=0)
noiseqData = readData(data = myTMM, factors=mydesign[,c(3,2)])
mycdTMM = dat(noiseqData, type = "cd", norm = T)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."
temp=dat2save(mycdTMM)
temp1=cbind(temp$DiagnosticTest[match(mydesign[,3],rownames(temp$DiagnosticTest)),3],mydesign[,2])
table(as.data.frame(temp1))
#      tumor normal  
#  FAILED 1  2
#  PASSED 97 90
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix,
	conditions=mydesign[mydesign[,3]%in%colnames(FilteredMatrix),c(3,2)]))
myDESEQ=counts(deseqFactors,normalized=T)
noiseqData = readData(data = myDESEQ, factors=mydesign[mydesign[,3]%in%colnames(FilteredMatrix),c(3,2)])
#mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
#explo.plot(mycountsbio, plottype = "boxplot",samples=NULL)#ya se nota la normalización
mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."
temp=dat2save(myDESEQ)
temp1=cbind(temp$DiagnosticTest[match(mydesign[,3],rownames(temp$DiagnosticTest)),3],mydesign[,2])
table(as.data.frame(temp1))
#      tumor normal  
#  FAILED 0  1
#  PASSED 98 91
myUQUA=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQUA, factors=mydesign[mydesign[,3]%in%colnames(FilteredMatrix),c(3,2)])
mycdUQUA = dat(noiseqData, type = "cd", norm = T)
temp=dat2save(mycdUQUA)
temp1=cbind(temp$DiagnosticTest[match(mydesign[,3],rownames(temp$DiagnosticTest)),3],mydesign[,2])
table(as.data.frame(temp1))
#      tumor normal  
#  FAILED 7  14
#  PASSED 91 78


myPCA = dat(noiseqData, type = "PCA", norm = TRUE, logtransf = TRUE)

par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "tissue.definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "tissue.definition")

