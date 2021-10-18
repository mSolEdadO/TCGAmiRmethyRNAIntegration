library(data.table)
library(caret)

files=list.files("/labs/csbig/multiomics/2021",full.names=T)
files=files[grep("eige",files)]
files=files[grep("tar",files,invert=T)]
data=lapply(files,fread)
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
data=do.call(cbind,data)
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
		function(x) t(data[x[1]:x[2],]))
	names(data)=c("CpGs","transcripts","miRNAs")
data=data$CpGs
results=nearZeroVar(data,saveMetrics= TRUE)
table(results$nzv)
#  FALSE
# 393132