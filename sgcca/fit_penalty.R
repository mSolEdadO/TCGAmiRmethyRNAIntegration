#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
penalty_cpgs=as.numeric(args[1])
penalty_transcripts=as.numeric(args[2])
penalty_mir=as.numeric(args[3])

library(mixOmics)

#take model descriptors<-----------------recycled
describe=function(data,pc,pt,pm){
	#subsample observations
	#data=lapply(data,function(x) x[sample(1:n,subn),])
	resus=wrapper.sgcca(data,penalty=c(pc,pt,pm),scale=T,
		scheme="centroid")
	#get results description
	description=as.data.frame(do.call(rbind,resus$AVE$AVE_X))
	description$nfeatures=sapply(resus$loadings,function(x) sum(x!=0))
	description$omic=rownames(description)
	description$penalty=resus$penalty
	colnames(description)[1]="AVE"
return(description)}

######################DATA TO LIST OF MATRIXES PER MOLECULAR LEVEL
library(data.table)
library(parallel)

files=list.files()
files=files[grep("eigeN",files)]
sizes=c(Basal=129,Her2=47,LumA=417,LumB=141,Normal=76)

cl <- parallel::makeCluster(10)
clusterExport(cl, c("describe","penalty_cpgs","penalty_transcripts",
	"penalty_mir","wrapper.sgcca","files","sizes","fread"))
resus=do.call(rbind,parLapply(cl,1:10,function(x) {
	#1 more than the samples coz of the rownames 
	i=lapply(sizes,function(x) c(1,sample(2:x,10)))
	data=lapply(1:5,function(x) fread(files[x],select=i[[x]]))
	data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
	data=do.call(cbind,data)
	#separate omics
	data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
		function(x) t(data[x[1]:x[2],]))
	names(data)=c("CpGs","transcripts","miRNAs")
	describe(data,penalty_cpgs,penalty_transcripts,penalty_mir)}))
stopCluster(cl)
file=paste(penalty_cpgs,penalty_transcripts,penalty_mir,"tsv",sep='.')
write.table(resus,file,sep='\t',quote=F,row.names=F)
