temp=read.table("penalty_search.tsv",sep='\t',header=T)

#######################################DIAGNOSTIC PLOTS
library(ggplot2)
library(gridExtra)

temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
temp=temp[order(temp$penalty),]
png("AVE.png",width=800)
 ggplot(temp,aes(y=AVE,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))
dev.off()
png("nfeatures.png",width=800)
 ggplot(temp,aes(y=nfeatures,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))+
 scale_y_continuous(trans="log10")
dev.off()
#plot median AVE vs meadian nfeatures
omics=levels(temp$omic)
omics=lapply(omics,function(x) temp[temp$omic==x,])
omics=lapply(omics,function(x) as.data.frame(apply(x,2,as.numeric)))
names(omics)=levels(temp$omic)
omics=lapply(omics,function(y) sapply(unique(y$penalty),function(x) 
apply(y[y$penalty==x,],2,median)))#better than mean?
omics=lapply(omics,function(x) as.data.frame(t(x)))
#indi plots or CpGs will determine axis
plots=lapply(1:3,function(x) ggplot(omics[[x]],
	aes(x=nfeatures,y=AVE,col=penalty))+geom_point()+
	ggtitle(names(omics)[x])+theme(text=element_text(size=18))+
	scale_x_continuous(trans="log10")+geom_line())
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

#######################################PENALTIES SUGGESTED BY PLOTS
grid=unique(temp$penalty)
#slopes=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
#	(y$AVE[x+1]-y$AVE[x])/(y$nfeatures[x+1]-y$nfeatures[x])))
slopes1=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
	abs(y$AVE[x+1]-y$AVE[x])/abs(y$nfeatures[x+1]-y$nfeatures[x])))

slopes1$miRNAs[abs(slopes1$miRNAs)=="Inf"]=NA
slopes1$transcripts[abs(slopes1$transcripts)=="Inf"]=NA
#penalty=sapply(slopes,function(x) grid[which.max(x)+1])
#       CpGs transcripts      miRNAs 
#       0.02        0.03        0.08 
#slopes1
#       CpGs transcripts      miRNAs 
#       0.02        0.10        0.08
#if u want the point before the biggest drop in AVE
#transcripts penalty should be 0.09 
penalty1[2]=0.09
#######################################FINAL SGCCA
library(igraph)
library(mixOmics)
library(data.table)

data=fread(paste(subtype,"normalized",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")

final=wrapper.sgcca(X=data,penalty=penalty1,scale=T,
	scheme="centroid",ncomp=20)#ncomp to explain 50% of transcripts matrix according to mfa.R
output=rbind(rowSums(do.call(rbind,final$AVE$AVE_X)),final$penalty)
rownames(output)=c("sum_AVE","penalty")

#####PLOT LOADINGS
temp=lapply(final$loadings,as.numeric)
temp=data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(names(temp)[x],temp[[x]]))))
colnames(temp)=c("omic","loading")
temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
png("loadings.png")
 ggplot(temp[temp$loadings!=0,],aes(x=omic,y=loadings))+
 geom_boxplot()+theme(text=element_text(size=18))
dev.off()

initial=wrapper.sgcca(X=final$X,penalty=rep(1,3),scale=T,
	scheme="centroid",ncomp=20)#ncomp to explain 50% of transcripts matrix according to mfa.R
rbind(rowSums(do.call(rbind,initial$AVE$AVE_X)),
	rowSums(do.call(rbind,final$AVE$AVE_X))) 
#          CpGs transcripts    miRNAs
#[1,] 0.7512414   0.5791811 0.5558149
#[2,] 0.6099801   0.5408933 0.5326386
temp$initial=unlist(lapply(initial$loadings,as.numeric))
plots=lapply(levels(temp$omic),function(x) ggplot(temp[temp$omic==x,],
	aes(y=final,x=initial))+geom_point()+ggtitle(x)+
    theme(text=element_text(size=18)))
png("loadings_change.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()