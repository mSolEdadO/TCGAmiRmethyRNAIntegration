load("subset.Rda")
library(mixOmics)
#(ggsci)#pa repetir la paleta de colores de MFA.R

#pimp data
data=list(mRNA=t(rna),miRNA=t(mirna),methylation=t(methy))
Y=design$subtype
design = matrix(1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
#1=correlation to maximize
diag(design) = 0


#Tuning the number of components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 60,design = design)
#colSums(do.call(cbind,sgccda.res$explained_variance))
#       mRNA       miRNA methylation           Y 
#  0.6959577   0.7630531   0.6887422  12.9489931 
#set.seed(123)# solo si no usas cpus
save(sgccda.res,file="splsda.Rda")
start_time <- Sys.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 5,cpus=5)#better with folds = 10 and  50-100 repeats, more time-consuming
end_time <- Sys.time()
end_time - start_time
#Time difference of 5.635757 hours
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]
#BER is appropriate in case of an unbalanced number of samples per class as it calculates the average proportion of wrongly classified samples in each class, weighted by the number of samples in each class
save(sgccda.res,perf.diablo,file="splsda.Rda")

test.keepX = list (mRNA = c(seq(10,100,10), 300, 500),
                    miRNA =  seq(10,140,10),
                    methylation = c(seq(10,100,10), seq(300, 1000, 200)))
#modify tune.block.splsda so parallel gets exported to every node in cl
#tmpfun <- get("tune.block.splsda", envir = asNamespace("mixOmics"))
#environment(tune.block.splsda) <- environment(tmpfun)
#attributes(tune.block.splsda) <- attributes(tmpfun)  # don't know if this is really needed
#assignInNamespace("tune.block.splsda", tune.block.splsda, ns="mixOmics")
start_time <- Sys.time()
seleccionadas = tune.block.splsda(X = sgccda.res$X, Y=sgccda.res$Y, ncomp = ncomp,test.keepX = test.keepX, 
                design = sgccda.res$design,validation = 'Mfold', folds = 10, nrepeat = 1,cpus = 4,
#for a more thorough tuning process, increase the nrepeat
                dist = "centroids.dist")
#This results in 2352 models being fitted for each component and each nrepeat
end_time <- Sys.time()
end_time - start_time
#Time difference of 1.148341 days
list.keepX = seleccionadas$choice.keepX

#final DIABLO model
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
#Time difference of 1.119338 mins
#

pdf("mixOmics.pdf")
#barplot(t(do.call(cbind,sgccda.final$explained_variance[1:3])*100),beside=T,las=2,col=scales::hue_pal()(4)[4:2],border=NA,ylab="% de varianza explicada",ylim=c(0,8))
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=1,title="metilación en comp1")
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=2,title="expresión de mRNA in comp1")
#plotLoadings(sgccda.final,comp=1,contrib="max",method="median",block=3,title="expresión de miRNA in comp1") 
circosPlot(sgccda.final, cutoff = 0.8, line =F,color.blocks=scales::hue_pal()(4)[4:2],comp=2,color.cor=c("chocolate3","grey20"))
plotArrow(sgccda.final,ind.names=F,col=ggsci::pal_lancet()(5)[sgccda.final$Y],legend=T,plot.arrows=F,pch=19,group=sgccda.final$Y)
dev.off()
