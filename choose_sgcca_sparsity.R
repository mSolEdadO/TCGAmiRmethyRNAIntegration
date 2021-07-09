#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(tidyverse)

######################NARROW DOWN SPARSITY VALUES
grid=seq(0,1,.2)
results=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=rep(x,3),scale=T))
#take model descriptors
AVE=as.data.frame(sapply(results,function(x) unlist(x$AVE)[1:3]))
nfeat=as.data.frame(sapply(results,function(x) 
	sapply(x$loadings,function(y) sum(y!=0))))
#format to plot
colnames(AVE)=grid
AVE$omic=rownames(nfeat)
AVE=AVE%>%pivot_longer(-7,names_to="sparsity",values_to="AVE")
colnames(nfeat)=grid
nfeat$omic=rownames(nfeat)
nfeat=nfeat%>%pivot_longer(-7,names_to="sparsity",
	values_to="nfeatures")
toplot=merge(nfeat,AVE,by=c("omic","sparsity"))
#force order
toplot$omic=factor(toplot$omic,levels=c("CpGs","transcripts","miRNAs"))
#plot
plots=lapply(levels(toplot$omic),function(i) 
	ggplot(toplot[toplot$omic==i,],aes(x=nfeatures,y=AVE))+
	geom_point()+ggtitle(i)+theme(text=element_text(size=18)))
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
 #no facet_wrap coz u want indy axes
dev.off()
#plot suggest CpG sparsity 0.2-0.4
#			  transcript   0-0.6
#			  miRNAs       0-0.2
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
conserved
#             CpGs       transcripts            miRNAs 
#     "cg26975351" "ENSG00000169116"    "hsa-mir-4756" 
temp=as.data.frame(sapply(results,function(x) sapply(1:3,function(y) 
	x$loadings[[y]][rownames(x$loadings[[y]])==conserved[y]])))
colnames(temp)=grid
temp$feature=conserved
temp=temp%>%pivot_longer(-7,names_to="sparsity",values_to="loading")
png("conserved_features.png")
 ggplot(temp,aes(x=sparsity,y=loading,group=feature))+
 geom_point(aes(color=feature))+theme(text=element_text(size=18))
dev.off()