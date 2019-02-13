load("MFA.Rda")
library(mixOmics)

#pimp data
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
#set.seed(123) solo si no usas cpus
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 5,cpus=5)
# better with folds = 10 and  50-100 repeats, more time-consuming
save(data,Y,sgccda.res,perf.diablo,file="splsda.Rda")

#Tuning keepX
test.keepX = list (mRNA = c(seq(10,100,10), 300, 500),
                   miRNA =  c(seq(10,100,10), seq(300, 1000, 200)),
                   methylation = c(seq(10,100,10), seq(300, 1000, 200)))
start_time <- Sys.time()
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")
#for a more thorough tuning process, increase the nrepeat
end_time <- Sys.time()
end_time - start_time
list.keepX = tune.TCGA$choice.keepX

#final DIABLO model
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

#load("MFA.Rda")
library(mixOmics)

#data to needed structure
#methy=do.call(rbind,subti$Her2[,1:395806])
#rna=do.call(rbind,subti$Her2[,395807:409710])
#mirna=do.call(rbind,subti$Her2[,409711:411298])
#colnames(methy)=rownames(subti$Her2)
#colnames(mirna)=rownames(subti$Her2)
#colnames(rna)=rownames(subti$Her2)
#data=list(mRNA=t(rna),miRNA=t(mirna),methylation=t(methy))
#Y=subti$Her2[,411300]
#names(Y)=rownames(subti$Her2)
#rm(MFAomics,MFAher2,subti,methy,rna,mirna)
#gc()
#design = matrix(1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(dat$
#1=correlation to maximize
#diag(design) = 0

#Tuning the number of components
#sgccda.res = block.splsda(X = data, Y = Y, ncomp = 60,design = design)
#save(data,Y,sgccda.res,file="splsda.Rda")
