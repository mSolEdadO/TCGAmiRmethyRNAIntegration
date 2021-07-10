#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(tidyverse)
library(pbapply)

######################NARROW DOWN SPARSITY VALUES
grid=seq(0,1,.2)
results=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=rep(x,3),scale=T))
#take model descriptors
describe=function(resus,names){
	AVE=as.data.frame(sapply(resus,function(x) 
		do.call(rbind,x$AVE[[1]])))
	nfeat=as.data.frame(sapply(resus,function(x) 
		sapply(x$loadings,function(y) sum(y!=0))))
	#format to plot
	colnames(AVE)=names
	AVE$omic=rownames(nfeat)
	AVE=AVE%>%pivot_longer(-omic,names_to="sparsity",values_to="AVE")
	colnames(nfeat)=names
	nfeat$omic=rownames(nfeat)
	nfeat=nfeat%>%pivot_longer(-omic,names_to="sparsity",
		values_to="nfeatures")
	toplot=merge(nfeat,AVE,by=c("omic","sparsity"))
	#force order
	toplot$omic=factor(toplot$omic,levels=c("CpGs","transcripts","miRNAs"))
	toplot$sparsity=as.numeric(toplot$sparsity)
return(toplot)}
description=describe(results,grid)
	#plot
plots=lapply(levels(description$omic),function(i) 
	ggplot(description[description$omic==i,],aes(x=nfeatures,
		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
	theme(text=element_text(size=18)))
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
 #no facet_wrap coz u want indy axes
dev.off()
#plot suggest CpG sparsity 0.2-0.4
#			  transcript   0-0.6
#			  miRNAs       0-0.2
#maybe turn it analytical by choosing the point when the slop changes
sapply(6:2,function(x) (description$AVE[x]-description$AVE[x-1])/
	(description$nfeatures[x]-description$nfeatures[x-1]))
#[1] 4.246797e-08 3.581813e-08 5.962024e-08 3.049028e-07 1.846931e-06

###############3######CHECK LOADINGS WITH NO SPARSITY
loadings=as.data.frame(do.call(rbind,results[[6]]$loadings))
loadings$omic=substr(rownames(loadings),1,1)
loadings$omic=gsub("E","transcripts",gsub("h","miRNAs",
	gsub("c","CpGs",loadings$omic)))
loadings$omic=factor(loadings$omic,
	levels=c("CpGs","transcripts","miRNAs"))
png("loadings_l1.png")
  ggplot(loadings,aes(y=comp1,x=omic))+geom_boxplot()+
  ylab("loadings")+theme(text=element_text(size=18))
dev.off()
#transcript & miRNA loadings are normally distributed
shapiro.test(loadings$comp1[loadings$omic=="transcripts"],5000))
#W = 0.99956, p-value = 0.3302
shapiro.test(loadings$comp1[loadings$omic=="miRNAs"])
#W = 0.99848, p-value = 0.8877
shapiro.test(sample(loadings$comp1[loadings$omic=="CpGs"],5000))#size must be between 3 and 5000
#W = 0.98626, p-value < 2.2e-16
#at sparsity=0, 1 feature is selected per omic
conserved=sapply(results[[1]]$loadings,function(x) rownames(x)[x!=0])
temp=as.data.frame(sapply(results,function(x) sapply(1:3,function(y) 
	x$loadings[[y]][rownames(x$loadings[[y]])==conserved[y]])))
colnames(temp)=grid
temp$feature=conserved
temp=temp%>%pivot_longer(-7,names_to="sparsity",values_to="loading")
png("conserved_features.png")
 ggplot(temp,aes(x=sparsity,y=loading,group=feature))+
 geom_point(aes(color=feature))+theme(text=element_text(size=18))
dev.off()

seq(0.05,0.7,length=6)
grid=seq(0.05,0.7,0.075)
results1=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(0.2,x,x),scale=T))
description=describe(results1,grid)
sapply(27:20,function(x) (description$AVE[x]-description$AVE[x-1])/(description$nfeatures[x]-description$nfeatures[x-1]))
#[1] 8.882339e-07 6.234892e-07 4.755818e-07 4.118815e-07 3.965331e-07
#[6] 4.858871e-07 9.089003e-07 3.766163e-06
#0.05 for transcripts
sapply(18:11,function(x) (description$AVE[x]-description$AVE[x-1])/(description$nfeatures[x]-description$nfeatures[x-1]))
#[1] -3.586748e-06 -2.357042e-06 -1.262106e-06 -1.070951e-06  1.677122e-06
#[6]  1.162186e-05  1.760406e-05  4.030400e-04
#0.2 for miRNAs
plots=lapply(levels(description$omic),function(i) 
	ggplot(description[description$omic==i,],aes(x=nfeatures,
		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
	theme(text=element_text(size=18)))
#png("sparsity_search1.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()
#plot suggest transcript sparsity 0.125-0.275
#			  miRNAs       		  0.2-0.275

final=wrapper.sgcca(data,penalty=c(0.2,0.05,0.2),scale=T))
