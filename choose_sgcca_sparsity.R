#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)
library(gridExtra)

#load sgcca.R output for a grid of values
files=list.files()
params=lapply(files,read.table,header=T)
grid=seq(.01,.9,length=10)#this was the grid used for sgcca.R

temp=do.call(rbind,params)
omics=levels(temp$omic)
sapply(omics,function(x) 
	cor.test(temp$explained_variance[temp$omic==x],
		temp$sparsity[temp$omic==x],method="spearman")$estimate)
#       CpGs.rho      miRNAs.rho transcripts.rho 
#      0.7537963       0.2520860       0.3010846 
#each has p-value < 2.2e-16 
#evar can only guide sparsity selection for CpGs
# cor.test(temp$explained_variance[temp$omic=="transcripts"],temp$sparsity[temp$omic=="CpGs"],method="spearman")
#S = 1.9659e+12, p-value < 2.2e-16

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

#sparsity of an omic affects selection of the other omics?
#separate per miR sparsity
pre_transcri=lapply(paste(grid,"params",sep='_'),function(x)
	do.call(rbind,pre[grep(x,names(pre))]))
#separate per transcript sparsity
pre_miR=lapply(paste("Her2",grid[5],grid,sep='_'),function(x)
	do.call(rbind,pre[grep(x,names(pre))]))
#test the effect of miR sparsity on transcript selection
sapply(2:10,function(x) wilcox.test(pre_transcri[[1]]$
	nfeatures[pre_transcri[[1]]$omic=="transcripts"],
	pre_transcri[[x]]$nfeatures[pre_transcri[[x]]$
	omic=="transcripts"],paired=T)$p.val)
#[1] 0.347945312 0.918818860 0.830158135 0.300657208 0.992860002
#[6] 0.664240708 0.308104697 0.009011566 0.362862150
#miR sparsity has no effect on transcript selection
sapply(2:10,function(x) wilcox.test(pre_miR[[1]]$
	nfeatures[pre_miR[[1]]$omic=="transcripts"],
	pre_miR[[x]]$nfeatures[pre_miR[[x]]$
	omic=="transcripts"],paired=T)$p.val)
#[1] 1.426783e-34 1.443684e-34 1.445942e-34 1.446452e-34 1.446655e-34
#[6] 1.446635e-34 1.446655e-34 1.007779e-34 3.719484e-36
#but transcript sparsity DOES affect miRNA selection

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


#será que la omica más grande determina a la que le sigue???
#cor(cpg sparsity, transcript evar)
lapply(omics[2:3],function(x) cor.test(temp$explained_variance[temp$omic==x],temp$sparsity[temp$omic==x],method="spearman"))