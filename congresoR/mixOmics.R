load("MFA.Rda")
library(mixOmics)

#data to needed structure
methy=do.call(rbind,subti$Her2[,1:395806])
rna=do.call(rbind,subti$Her2[,395807:409710])
mirna=do.call(rbind,subti$Her2[,409711:411298])
colnames(methy)=rownames(subti$Her2)
colnames(mirna)=rownames(subti$Her2)
colnames(rna)=rownames(subti$Her2)
data=list(mRNA=t(rna),miRNA=t(mirna),methylation=t(methy))
Y=subti$Her2[,411300]
names(Y)=rownames(subti$Her2)
rm(MFAomics,MFAher2,subti,methy,rna,mirna)
gc()
design = matrix(1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
#1=correlation to maximize
diag(design) = 0


#Tuning the number of components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 60,design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 100,cpus=7)
# better with folds = 10 and higher nrepeat value, but more time-consuming
save(data,Y,sgccda.res,perf.diablo,file="splsda.Rda")