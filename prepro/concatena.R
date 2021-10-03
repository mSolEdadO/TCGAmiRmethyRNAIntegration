#!/usr/bin/env Rscript

library(data.table)

subtype=read.table("subtype.tsv",header=T,sep='\t')
#revert comments for unnormalized data
#expre=fread("RNAseqnormalized.tsv")
#miR=fread("/miRNAseqNormi.tsv")
#methy=fread("methyM.tsv")
files=list.files("/labs/csbig/multiomics/2021",full.names=T)
files=files[grep("eigenNormi",files)]
data=lapply(files,fread)
#expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
#miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
#methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(data)=gsub(".eigenNormi","",files)
names(data)=gsub("RNAseqnormalized","transcripts",
	gsub("miRNAseqNormi","miRNAs",gsub("methyM","CpGs",names(data))))
print(sapply(data,dim))
#print(sapply(data,function(x) head(rownames(x))))

#choose methy order
#subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#expre=expre[,order(match(colnames(expre),subtype$samples))]
#miR=miR[,order(match(colnames(miR),subtype$samples))]
subtype=subtype[order(match(subtype$samples,colnames(data$CpGs))),]
data[2:3]=lapply(data[2:3],function(x)
	x[,order(match(colnames(x),subtype$samples))])
print(names(data))

#data per subtype
concatenated=lapply(levels(subtype$subtype),function(x) 
	list(CpGs=data$CpGs[,subtype$subtype==x],
		transcripts=data$transcripts[,subtype$subtype==x],
		miRNA=data$miRNAs[,subtype$subtype==x]))
names(concatenated)=levels(subtype$subtype)

##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
print(sapply(concatenated,dim))
#      Basal   Her2   LumA   LumB Normal
#[1,] 410813 410813 410813 410813 410813
#[2,]    128     46    416    140     75
lapply(1:5,function(x) write.table(concatenated[[x]],
	paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))

#to go back
#apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
#	function(x) data[x[1]:x[2],])

#next check PCs per subtype & data→→→→→→→→→→→→→→mfa.R
#next check inter-omic correlation→→→→→→→corrKernel.R<----optional
#drop near zero var features???????

###########################
#join omics
#!/usr/bin/env Rscript
#library(data.table)

#files=list.files()
#files=files[grep("mtrx",files)]
#data=lapply(files,fread)
#data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
#data=do.call(cbind,data)
#data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
# function(x) data[x[1]:x[2],])
#names(data)=c("CpGs","transcripts","miRNAs")
#write.table(data[[1]],paste(names(data)[1],"mtrx",sep='.'),
#		sep='\t',quote=F)

