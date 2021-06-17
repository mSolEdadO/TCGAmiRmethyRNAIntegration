#!/usr/bin/env Rscript

library(data.table)
library(mixKernel)
between_kernel_corr=function(data,subtype){
	#compute linear kernel, good performance with lots of features
	#nonlinear kernels are not necessarily more accurate
	#https://stats.stackexchange.com/questions/73032/linear-kernel-and-non-linear-kernel-for-support-vector-machine
	k=lapply(data,compute.kernel)
#no_cores <- detectCores() - 1  
#registerDoParallel(cores=no_cores)  
#cl <- makeCluster(no_cores, type="FORK") 
#k=parLapply(cl,data,compute.kernel)
	#plot between kernel corr
	png(paste("corr",subtype,"png",sep='.'))
	print({cim.kernel(CpGs=k$CpGs,transcripts=k$transcripts,
	 	miRNAs=k$miRNAs,method=c("shade"))})
	dev.off()
	cim.kernel(CpGs=k$CpGs,transcripts=k$transcripts,
		miRNAs=k$miRNAs,method=c("number"))
}

subtype = commandArgs(trailingOnly=TRUE)
#load data
data=fread(paste(subtype,"mtrx",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1))
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),
#numbers change with the number of features
	1,function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
print(subtype)
#analyze subtype
between_kernel_corr(data,subtype)


#####################################RESULTS
#[1] "Her2"
#                 CpGs transcripts    miRNAs
#CpGs        1.0000000   0.7219649 0.7089221
#transcripts 0.7219649   1.0000000 0.9065280
#miRNAs      0.7089221   0.9065280 1.0000000
#[1] "LumA"
#                 CpGs transcripts    miRNAs
#CpGs        1.0000000   0.2835642 0.2522872
#transcripts 0.2835642   1.0000000 0.5190197
#miRNAs      0.2522872   0.5190197 1.0000000
#[1] "LumB"
#                 CpGs transcripts    miRNAs
#CpGs        1.0000000   0.4781921 0.4434121
#transcripts 0.4781921   1.0000000 0.7703846
#miRNAs      0.4434121   0.7703846 1.0000000
#[1] "Basal"
#                 CpGs transcripts    miRNAs
#CpGs        1.0000000   0.5156423 0.4759017
#transcripts 0.5156423   1.0000000 0.7763678
#miRNAs      0.4759017   0.7763678 1.0000000
#[1] "Normal"
#                 CpGs transcripts    miRNAs
#CpGs        1.0000000   0.3617339 0.4723301
#transcripts 0.3617339   1.0000000 0.7019761
#miRNAs      0.4723301   0.7019761 1.0000000
