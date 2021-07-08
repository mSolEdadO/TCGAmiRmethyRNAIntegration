#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)
library(gridExtra)

#load sgcca.R output for a grid of values
files=list.files()
params=lapply(files,read.table,header=T)
grid=seq(.01,.9,length=10)#this was the grid used for sgcca.R

####################DESCRIPTION OF DATA##########################
temp=do.call(rbind,params)
omics=levels(temp$omic)
sapply(omics,function(x) 
	cor.test(temp$explained_variance[temp$omic==x],
		temp$sparsity[temp$omic==x],method="spearman")$estimate)
#       CpGs.rho      miRNAs.rho transcripts.rho 
#      0.7537963       0.2520860       0.3010846 
#each has p-value < 2.2e-16 
#evar can only guide sparsity selection for CpGs
#cor.test(temp$explained_variance[temp$omic=="transcripts"],temp$sparsity[temp$omic=="CpGs"],method="spearman")
#S = 1.9659e+12, p-value < 2.2e-16

########sparsity of an omic affects selection of the other omics?
#separate per miR sparsity
i=as.character(sapply(grid,function(x) 
	sapply(grid,function(y) paste("Her2",x,y,sep='_'))))
temp=lapply(i,function(x) params[grep(x,files)])
#effect of varying miR sparsity on transcript selection
table(unlist(lapply(temp,function(x) sapply(1:9,function(y) sapply((y+1):10,function(z) wilcox.test(x[[y]]$nfeatures[x[[y]]$omic=="transcripts"],x[[z]]$nfeatures[x[[z]]$omic=="transcripts"],paired=T)$p.val))))<0.05,useNA="always")
#FALSE  TRUE  <NA> 
# 3860   189   451 
#most comparisons are not significant, miRNA sparsity has little effect
#effect of varying miRNA sparsity on CpG selection
table(unlist(lapply(temp,function(x) sapply(1:9,function(y) sapply((y+1):10,function(z) wilcox.test(x[[y]]$nfeatures[x[[y]]$omic=="CpGs"],x[[z]]$nfeatures[x[[z]]$omic=="CpGs"],paired=T)$p.val))))<0.05,useNA="always")
#FALSE  TRUE  <NA> 
# 3572   478   450
#miR sparsity affects more CpG selection than transcript selection

#separate per CpG sparsity
i=as.character(sapply(grid,function(x) 
	sapply(grid,function(y) paste(x,y,"params",sep='_'))))
temp=lapply(i,function(x) params[grep(x,files)])
#effect of varying CpG sparsity on transcript selection
table(unlist(lapply(temp,function(x) sapply(1:9,function(y) sapply((y+1):10,function(z) wilcox.test(x[[y]]$nfeatures[x[[y]]$omic=="transcripts"],x[[z]]$nfeatures[x[[z]]$omic=="transcripts"],paired=T)$p.val))))<0.05,useNA="always")
#FALSE  TRUE  <NA> 
# 3521   529   450 
#CpG sparsity effect is larger than miR sparsity effect
#effect of varying CpG sparsity on miRNA selection
table(unlist(lapply(temp,function(x) sapply(1:9,function(y) sapply((y+1):10,function(z) wilcox.test(x[[y]]$nfeatures[x[[y]]$omic=="miRNAs"],x[[z]]$nfeatures[x[[z]]$omic=="miRNAs"],paired=T)$p.val))))<0.05,useNA="always")
#FALSE  TRUE  <NA> 
# 2867   729   904 
#CpG sparsity affects more miR selection than transcript selection
######################separate per transcript sparsity
i=as.character(sapply(grid,function(x) sapply(grid,function(y) paste("Her2",x,"0.*?",y,"params",sep='_'))))
temp=lapply(intersect(i,j),function(x) params[grep(x,files)])
#thus u CAN plot all together 

####################CHOOSING SPARSITY##########################
#this ignores the sparsity value of non-ploted omics
plotit=function(data,level){
	data$sparsity=as.character(data$sparsity)
	p=ggplot(data[data$omic==level,],
	aes(x=sparsity,y=explained_variance))+geom_boxplot()+
	scale_x_discrete(labels=signif(grid,digits=2))+
	theme(text=element_text(size=18))
	p1=ggplot(data[data$omic==level,],aes(x=sparsity,y=nfeatures))+
	geom_boxplot()+#scale_y_continuous(trans="log10")+
	scale_x_discrete(labels=signif(grid,digits=2))+
	theme(text=element_text(size=18))
return(plots=list(p,p1))}
plots=plotit(temp,"CpGs")
png("CpG_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()
#so far, 0.5 or →0.4← for CpG sparsity

#CpG sparsity fixed at 0.4
i=grep("Her2_0.4",files)
pre=params[i]
names(pre)=files[i]

temp=do.call(rbind,pre)
#cor.test(temp$explained_variance[temp$omic=="miRNAs"],temp$sparsity[temp$omic=="transcripts"],method="spearman")
#S = 1350124768, p-value = 0.5735
plots=plotit(temp,"transcripts")
png("transcript_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()
#so far, 0.41 or →0.31← for transcript sparsity

pre=pre[grep(paste("_",grid[4],"_0",sep=''),names(pre))]
temp=do.call(rbind,pre)
plots=plotit(temp,"miRNAs")
png("miRNA_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()
#so far, 0.7 or →0.6← for transcript sparsity

#check selection
temp=params[[grep(paste("Her2",grid[5],grid[4],grid[7],sep='_'),
	files)]]
sapply(omics,function(x) colMeans(temp[temp$omic==x,1:2]))
#                           CpGs       miRNAs  transcripts
#explained_variance 1.372375e-01   0.05925499 6.288162e-02
#nfeatures          9.947385e+04 374.25000000 2.878450e+03
################################
nopen=wrapper.sgcca(data,penalty=c(1,1,1),scale=T)
loadings=as.data.frame(do.call(rbind,nopen$loadings))
loadings$omic=substr(rownames(loadings),1,1)
loadings$omic=gsub("E","transcripts",gsub("h","miRNAs",
	gsub("c","CpGs",loadings$omic)))
loadings$omic=factor(loadings$omic,
	levels=c("CpGs","transcripts","miRNAs"))
png("loadings_l0.png")
  ggplot(loadings,aes(y=comp1,x=omic))+geom_boxplot()+
  ylab("loadings")+theme(text=element_text(size=18))
dev.off()
#transcript & miRNA loadings are normal distributed
shapiro.test(sample(nopen$loadings$transcripts,5000))
#W = 0.99956, p-value = 0.3302
shapiro.test(nopen$loadings$miRNAs)
#W = 0.99848, p-value = 0.8877
shapiro.test(sample(nopen$loadings$CpGs,5000))
#W = 0.98626, p-value < 2.2e-16
