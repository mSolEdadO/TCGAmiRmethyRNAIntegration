#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(tidyverse)
library(pbapply)

######################NARROW DOWN SPARSITY VALUES
grid=seq(0,.9,.1)
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(x,1,1),scale=T))
#take model descriptors<-----------------recycled
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

results=list()
results$CpGs=describe(temp,grid)
results$CpGs=results$CpGs[results$CpGs$omic=="CpGs",]
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,x,1),scale=T))
results$transcripts=describe(temp,grid)
results$transcripts=results$transcripts[results$
	transcripts$omic=="transcripts",]
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,1,x),scale=T))
results$miRNAs=describe(temp,grid)
results$miRNAs=results$miRNAs[results$
	miRNAs$omic=="miRNAs",]
results=do.call(rbind,results)

plots=lapply(levels(results$omic),function(i) 
	ggplot(results[results$omic==i,],aes(x=nfeatures,
		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
	theme(text=element_text(size=18)))
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
 #no facet_wrap coz u want indy axes
dev.off()
#turn it analytical by choosing the point when the slop changes????
sapply(10:2,function(x) (description$AVE[x]-description$AVE[x-1])/
	(description$nfeatures[x]-description$nfeatures[x-1]))
#[1] 4.265165e-08 3.532026e-08 3.315587e-08 4.461965e-08 5.877692e-08
#[6] 7.409449e-08 2.620996e-07 4.787702e-07 9.293855e-06
#sparsity: 0.2-0.3 for CpGs
sapply(20:12,function(x) (results$AVE[x]-results$AVE[x-1])/(results$nfeatures[x]-results$nfeatures[x-1]))
#[1]  7.945155e-08  8.699518e-08  1.455164e-07  2.795368e-07  7.898728e-07
#[6]  1.674634e-06  1.457345e-06 -2.883192e-06  5.426026e-06
#sparsity: 0.6-0.7 for transcripts
sapply(30:22,function(x) (results$AVE[x]-results$AVE[x-1])/(results$nfeatures[x]-results$nfeatures[x-1]))
#[1] -9.330695e-08  1.483156e-07  4.117869e-07  3.723666e-06  5.615809e-06
#[6]  1.169444e-05  2.834565e-05  4.721406e-05  2.599923e-04
#sparsity: 0.5-0.6 for miRNAs

#####################CHECK LOADINGS
init=wrapper.sgcca(data,penalty=c(1,1,1),scale=T)
final=wrapper.sgcca(data,penalty=c(0.2,0.6,0.5),scale=T)

loadings=as.data.frame(do.call(rbind,init$loadings))
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
shapiro.test(sample(loadings$comp1[loadings$omic=="transcripts"],5000))
#W = 0.99949, p-value = 0.1937
shapiro.test(loadings$comp1[loadings$omic=="miRNAs"])
#W = 0.99848, p-value = 0.8877
shapiro.test(sample(loadings$comp1[loadings$omic=="CpGs"],5000))#size must be between 3 and 5000
#W = 0.98626, p-value < 2.2e-16
###final loadgings all have p-value < 2.2e-16

################################COMPARE INIT & FINAL
#%of selected features
100*sapply(final$loadings,function(x) 
	sum(x!=0))/sapply(init$loadings,function(x) sum(x!=0))
#       CpGs transcripts      miRNAs 
#   6.187235   63.565029   46.026490 
#% of AVE from initial
100*do.call(rbind,final$AVE[[1]])/do.call(rbind,init$AVE[[1]])
#                comp1
#CpGs         80.10103
#transcripts  96.51604
#miRNAs      103.08098
#####whats with miRNAs?????
temp=as.data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(init$loadings[[x]],final$loadings[[x]]))))
colnames(temp)=c("initial","final")
temp$omic=substr(rownames(temp),1,1)
temp$omic=gsub("E","transcripts",gsub("h","miRNAs",
 gsub("c","CpGs",loadings$omic)))
plots=lapply(unique(temp$omic),function(x) 
	ggplot(temp[temp$omic==x,],aes(x=initial,y=final))+geom_point()+
	theme(text=element_text(size=18)))
png("loadings_change.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()