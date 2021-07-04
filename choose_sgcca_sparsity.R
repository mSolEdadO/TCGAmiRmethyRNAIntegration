#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
library(ggplot2)

#load sgcca.R output for a grid of values
files=list.files()
params=lapply(files,read.table,header=T)

temp=do.call(rbind,params)
omics=levels(temp$omic)
sapply(omics,function(x) 
	cor.test(temp$explained_variance[temp$omic==x],
		temp$sparsity[temp$omic==x],method="spearman")$estimate)
#       CpGs.rho      miRNAs.rho transcripts.rho 
#      0.7537963       0.2520860       0.3010846 
#evar can only guide sparsity selection for CpGs

#this ignores the sparsity value of non-ploted omics
plotit=function(data,level){
	data$sparsity=as.character(data$sparsity)
	p=ggplot(temp[temp$omic==level,],
	aes(x=sparsity,y=explained_variance))+geom_boxplot()+
	scale_x_discrete(labels=signif(as.numeric(unique(temp$sparsity)),
	digits=2))+theme(text=element_text(size=18))
	p1=ggplot(temp[temp$omic==level,],aes(x=sparsity,y=nfeatures))+
	geom_boxplot()+scale_y_continuous(trans="log10")+
	scale_x_discrete(labels=signif(as.numeric(unique(temp$sparsity)),
	digits=2))+theme(text=element_text(size=18))
return(plots=list(p,p1))}
plots=plotit(temp,"CpGs")
png("CpG_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()
#so far, 0.5 or →0.4← for CpG sparsity

#CpG sparsity fixed at 0.4
grid=seq(.01,.9,length=10)#this was the grid used for sgcca.R
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
plots=plotit(temp,"transcripts")
png("transcript_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()
#so far, 0.31 or →0.21← for transcript sparsity

pre=pre[grep(paste("_",grid[3],"_0",sep=''),names(pre))]
temp=do.call(rbind,pre)
plots=plotit(temp,"miRNAs")
png("miRNA_sparsity.png")
 grid.arrange(plots[[1]],plots[[2]])
dev.off()


#será que la omica más grande determina a la que le sigue???
#cor(cpg sparsity, transcript evar)