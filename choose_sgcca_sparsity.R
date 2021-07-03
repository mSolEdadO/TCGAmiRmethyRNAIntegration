#load sgcca.R output for a grid of values
files=list.files()
params=lapply(files,read.table,header=T)

#sparsity of the other omics matters?
temp=do.call(rbind,params[grep("Her2_0.01_0.01",files)])
#CpG s == transcript s == 0.01, all miRNA s values
temp1=do.call(rbind,params[grep("Her2_0.9_0.9",files)])
#CpG s == transcript s == 0.9, all miRNA s values
wilcox.test(temp$nfeatures[temp$omic=="miRNAs"],
	temp1$nfeatures[temp1$omic=="miRNAs"],paired=T)
#V = 8385.5, p-value = 1.336e-13
#it does: distribution of selected miRNAs changes between
# the 2 sparsities of CpGs & transcripts 
