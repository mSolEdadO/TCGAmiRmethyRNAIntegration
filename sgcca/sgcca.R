#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
ncomp=args[2]

library(igraph)
library(mixOmics)
library(data.table)
########################DATA
data=fread(paste("/labs/csbig/multiomics/2021",
	paste(subtype,"eigeNormi",sep='.'),sep='/'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
penalty=c(CpGs=0.02,transcripts=0.02,miRNAs=0.05)#output of choose_penalty.R

########################THE SGCCA
final=wrapper.sgcca(X=data,penalty=penalty,scale=F,
	scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
save(final,file=paste(subtype,"RData",sep='.'))

#####PLOT LOADINGS
temp=lapply(final$loadings,as.numeric)
temp=data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(names(temp)[x],temp[[x]]))))
colnames(temp)=c("omic","loading")
temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
png(paste(subtype,"loadings.png",sep='-'))
 ggplot(temp[temp$loadings!=0,],aes(x=omic,y=loadings))+
 geom_boxplot()+theme(text=element_text(size=18))
dev.off()

initial=wrapper.sgcca(X=data,penalty=rep(1,3),scale=F,
	scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
rbind(rowSums(do.call(rbind,initial$AVE$AVE_X)),
	rowSums(do.call(rbind,final$AVE$AVE_X))) 
#          CpGs transcripts    miRNAs
#[1,] 0.7512414   0.5791811 0.5558149
#[2,] 0.6099801   0.5408933 0.5326386
temp$initial=unlist(lapply(initial$loadings,as.numeric))
plots=lapply(levels(temp$omic),function(x) ggplot(temp[temp$omic==x,],
	aes(y=final,x=initial))+geom_point()+ggtitle(x)+
    theme(text=element_text(size=18)))
png(paste(subtype,"loadings_change.png",sep='-'))
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

