get.minRMSE=function(folder){
	files=list.files(folder,full.names=T)
	files=files[grep("results",files)]
	files=lapply(unique(sapply(strsplit(files,".",fixed=T),
		function(x) x[4])),function(y) files[grep(y,files)])
	modelos=lapply(files,function(x) lapply(x,function(y) read.table(y,sep='\t',header=T)))
	names(modelos[[1]])=sapply(strsplit(sapply(strsplit(files[[1]],"/"),function(x) 
		x[2]),'.',fixed=T),function(y) y[1])
	names(modelos[[2]])=names(modelos[[1]])
	names(modelos[[3]])=names(modelos[[1]])
	names(modelos[[4]])=names(modelos[[1]])
	names(modelos[[5]])=names(modelos[[1]])
	names(modelos)=sapply(1:5,function(z) 
		sapply(strsplit(sapply(strsplit(files[[z]][1],"/"),function(x)
			x[2]),".",fixed=T),function(y) y[2]))
	return(lapply(modelos,function(x) sapply(x,function(y) y[which(y$RMSE==min(y$RMSE))[1],])))
	}

#get characteristics of the best model per scheme and gene
cg01ENSG1hsa05=get.minRMSE("cg0.1ENSG1hsa0.5")
cg01ENSG1hsa01=get.minRMSE("cg0.1ENSG1hsa0.1")
cg05ENSG1hsa01=get.minRMSE("cg0.5ENSG1hsa0.1")
cg1ENSG1hsa1=get.minRMSE("cg1.0ENSG1hsa1.0")
opci=list(cg1ENSG1hsa1,cg01ENSG1hsa01,cg01ENSG1hsa05,cg05ENSG1hsa01)
names(opci)=list("cg1ENSG1hsa1","cg01ENSG1hsa01","cg01ENSG1hsa05","cg05ENSG1hsa01")

sapply(opci,function(x) sapply(x,function(y) median(as.numeric(y[3,]),na.rm=T)))
lapply(1:5,function(z) table(apply(do.call(rbind,lapply(1:4,function(x) opci[[x]][[z]][3,])),2,function(y) names(opci)[y==min(as.numeric(y))])))
temp=lapply(1:50,function(z) unlist(lapply(1:4,function(x) median(sapply(1:5,function(y) as.numeric(opci[[x]][[y]][3,z]))))))
lapply(temp,function(x) names(opci)[x==min(x)])