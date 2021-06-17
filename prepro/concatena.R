library(data.table)

subtype=read.table("subtype.tsv",header=T,sep='\t')
expre=fread("parallel-aracne/RNAseqnormalized.tsv")
miR=fread("parallel-aracne/miRNAseqNormi.tsv")
methy=fread("methyM.tsv")
expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
#choose methy order
subtype=subtype[order(match(subtype$samples,colnames(methy))),]
expre=expre[,order(match(colnames(expre),subtype$samples))]
miR=miR[,order(match(colnames(miR),subtype$samples))]

#data per subtype
concatenated=lapply(levels(subtype$subtype),function(x) 
	list(methy=methy[,subtype$subtype==x],
		RNA=expre[,subtype$subtype==x],
		miRNA=miR[,subtype$subtype==x]))
names(concatenated)=levels(subtype$subtype)

##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
sapply(concatenated,dim)
#      Basal   Her2   LumA   LumB Normal
#[1,] 410813 410813 410813 410813 410813
#[2,]    128     46    416    140     75
lapply(1:5,function(x) write.table(concatenated[[x]],
	paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))

#to go back
#apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
#	function(x) data[x[1]:x[2],])

#next check PCs per subtype & data→→→→→→→→→→→→→→mfa.R
#next check inter-omic correlation→→→→→→→corrKernel.R
#drop near zero var features???????''

