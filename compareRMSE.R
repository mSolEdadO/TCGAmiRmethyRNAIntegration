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
cg05ENSG1hsa01=get.minRMSE("cg0.5ENSG1hsa0.1")
cg01ENSG1hsa01=get.minRMSE("cg0.1ENSG1hsa0.1")
cg05ENSG1hsa05=get.minRMSE("cg0.5ENSG1hsa0.5")
cg1ENSG1hsa1=get.minRMSE("cg1.0ENSG1hsa1.0")
opci=list(cg1ENSG1hsa1,cg01ENSG1hsa01,cg01ENSG1hsa05,cg05ENSG1hsa01,cg05ENSG1hsa05)
names(opci)=list("cg1ENSG1hsa1","cg01ENSG1hsa01","cg01ENSG1hsa05","cg05ENSG1hsa01","cg05ENSG1hsa05")

elegidos=do.call(cbind,lapply(1:5,function(x)
 apply(sapply(opci,function(y) y[[x]][3,]),1,function(z) names(opci)[z==min(as.numeric(z))])))
colnames(elegidos)=names(opci$cg1ENSG1hsa1)
lambdas=t(sapply(1:50,function(x) sapply(1:5,function(y) opci[[elegidos[x,y]]][[y]][,x]$lambda)))
modelos=as.character(sapply(1:50,function(x) sapply(1:5,function(y)
 paste(colnames(elegidos)[y],rownames(elegidos)[x],lambdas[x,y],elegidos[x,y],sep=' '))))
modelos=gsub("hsa"," ",gsub("cg","",gsub("hsa0"," 0.",gsub("ENSG1"," 1",gsub("cg0","0.",modelos)))))
writeLines(modelos,"bestModels.tsv")