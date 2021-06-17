#!/usr/bin/env Rscript

between_kernel_corr=function(data,subtype){
	print(subtype)
	#compute linear kernel
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
	p=cim.kernel(CpGs=k$CpGs,transcripts=k$transcripts,
		miRNAs=k$miRNAs,method=c("number"))
	print("Correlation between kernels")
	print(p)
}

subtype = commandArgs(trailingOnly=TRUE)
#load data
data=fread(paste(subtype,"mtrx",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1))
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),
#numbers change with the number of features
	1,function(x) t(data[[2]][x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
print(subtype)
#analyze subtype
between_kernel_corr(data,subtype)
