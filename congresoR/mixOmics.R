load("MFA.Rda")
library(mixOmics)
library(ggsci)#pa repetir la paleta de colores de MFA.R

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
#perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 5,cpus=5)
# better with folds = 10 and  50-100 repeats, more time-consuming
#save(data,Y,sgccda.res,perf.diablo,file="splsda.Rda")

#Tuning keepX
#test.keepX = list (mRNA = c(seq(10,100,10), 300, 500),
#                   miRNA =  c(seq(10,100,10), seq(300, 1000, 200)),
#                   methylation = c(seq(10,100,10), seq(300, 1000, 200)))
#start_time <- Sys.time()
#tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
#                              test.keepX = test.keepX, design = design,
#                              validation = 'Mfold', folds = 10, nrepeat = 1,
#                              cpus = 2, dist = "centroids.dist")
#for a more thorough tuning process, increase the nrepeat
#end_time <- Sys.time()
#end_time - start_time
#list.keepX = tune.TCGA$choice.keepX

#final DIABLO model
#sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
#                          keepX = list.keepX, design = design)
#
#
sgccda.res1 = block.splsda(X = data, Y = Y, ncomp =55,
	keepX=list(mRNA = rep(40,55), miRNA = rep(20,55),methylation = rep(20,55)),
	design=sgccda.res$design,max.iter=1000)

pdf("mixOmics.pdf")
#barplot(t(do.call(cbind,sgccda.final$explained_variance[1:3])*100),beside=T,las=2,col=scales::hue_pal()(4)[4:2],border=NA,ylab="% de varianza explicada",ylim=c(0,8))
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=1,title="metilación en comp1")
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=2,title="expresión de mRNA in comp1")
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=3,title="expresión de miRNA in comp1") 
circosPlot(sgccda.final, cutoff = 0.8, line =F,color.blocks=scales::hue_pal()(4)[4:2],comp=2,color.cor=c("chocolate3","grey20"))
plotArrow(sgccda.final,ind.names=F,col=ggsci::pal_lancet()(5)[sgccda.final$Y],legend=T,plot.arrows=F,pch=19,group=sgccda.final$Y)
dev.off()

selected=sapply(1:55,function(y) unlist(sapply(selectVar(sgccda.res1,comp=y)[1:3],function(x) x[[1]])))
selected=lapply(c("ENSG","hsa","cg"),function(x) selected[grep(x,selected)])
subset=lapply(1:3,function(x) data[[x]][,colnames(data[[x]])%in%selected[[x]]])
names(subset)=c("rna","mirna","methy")
library(r.jive)
JIVEsubset=jive(lapply(subset,function(x) t(x)))
library(FactoMineR)
MFAsubset=MFA(scale(do.call(cbind,subset)),group=c(2132,895,1100),name.group=c("rna","mirna","methy"),ncp=5,graph=F)



################subset237###############
start_time <- Sys.time()
tuneando = tune.block.splsda(X = data, Y = Y, ncomp = 5,test.keepX = test.keepX, design = design,validation = 'Mfold', folds = 10, nrepeat = 1,cpus = 4, dist = "centroids.dist")
end_time <- Sys.time()
end_time - start_time
#Time difference of 27.37144 mins

pdf("mixOmics237.pdf")
barplot(t(do.call(cbind,sgccda.res$explained_variance[1:3])*100),beside=T,las=2,col=scales::hue_pal()(4)[4:2],border=NA,ylab="% of omics variance explained",ylim=c(0,70))
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE,col=scales::hue_pal()(4)[2:4],pch=rep(19,3))
plotLoadings(sgccda.res,contrib="max",method="median",block=1,title="mRNA in comp1",legend.color=ggsci::pal_lancet()(5),comp=1)
plotLoadings(sgccda.res,contrib="max",method="median",block=2,title="miRNA in comp3",legend.color=ggsci::pal_lancet()(5),comp=3)
plotLoadings(sgccda.res,contrib="max",method="median",block=3,title="methylation in comp2",legend.color=ggsci::pal_lancet()(5),comp=2)
circosPlot(sgccda.res, cutoff = 0.7, line = T,size.labels=2,color.blocks=scales::hue_pal()(4)[4:2],,color.cor=c("chocolate3","grey20"),color.Y=ggsci::pal_lancet()(5))
dev.off()
