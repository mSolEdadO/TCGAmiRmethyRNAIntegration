library(ggplot2)

#load sgcca.R output for a grid of values
files=list.files()
params=lapply(files,read.table,header=T)

temp=do.call(rbind,params)
temp$sparsity=as.character(temp$sparsity)
temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
png("evar.png",width=1000)
	ggplot(temp,aes(x=sparsity,y=explained_variance))+
	geom_boxplot()+facet_wrap(~omic)+
	scale_x_discrete(labels=signif(as.numeric(unique(temp$sparsity)),
	digits=2))+theme(text=element_text(size=18))
dev.off()
#evar seems to be indi of sparsity, for miR & transcri
png("nfeatures.png",width=1000)
	ggplot(temp,aes(x=sparsity,y=nfeatures))+
	geom_boxplot()+scale_y_continuous(trans="log10")+facet_wrap(~omic)+
	scale_x_discrete(labels=signif(as.numeric(unique(temp$sparsity)),
		digits=2))+theme(text=element_text(size=18))
dev.off()
#lower s, wider nfeatures but with smaller values
#ideally max evar & min nfeatures

#sparsity of the other omics matters?
temp=do.call(rbind,params[grep("Her2_0.01_0.01",files)])
#CpG s == transcript s == 0.01, all miRNA s values
temp1=do.call(rbind,params[grep("Her2_0.9_0.9",files)])
#CpG s == transcript s == 0.9, all miRNA s values
wilcox.test(temp$nfeatures[temp$omic=="miRNAs"],
	temp1$nfeatures[temp1$omic=="miRNAs"],paired=T)
#V = 8385.5, p-value = 1.336e-13
#it does: distri of selected miRs changes between the 2 sparsities
#this may explain evar.png