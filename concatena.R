load("ini/subtiTMMArsyn.RData")
load("ini/mirTMMARSyn.RData")
load("ini/imptMethy.RData")

shared=lapply(c("LumA","Basal","LumB","Her2","normal"),function(x) 
			sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
				substr(colnames(y[[which(names(y)==x)]]),1,12)))
shared=sapply(1:5,function(x) 
			intersect(intersect(shared[[x]][[1]],shared[[x]][[2]]),shared[[x]][[3]]))
names(shared)=c("LumA","Basal","LumB","Her2","normal")
concatenadas=lapply(1:5,function(x) sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
	y[[which(names(y)==names(shared)[x])]][,substr(colnames(y[[which(names(y)==names(shared)[x])]]),1,12)%in%shared[[x]]]))
names(concatenadas)=names(shared)
concatenadas=lapply(concatenadas,function(x) sapply(x,function(y) y[,order(substr(colnames(y),1,12))]))
sapply(concatenadas,function(x) ncol(x[[1]]))
#  LumA  Basal   LumB   Her2 normal 
#   331    135    177     75     75 

save(concatenadas,methyDesign,subtipos,file="porSubti.RData")
#########################################
#lapply(concatenadas,function(x) sapply(x,function(y) summary(as.numeric(y))))
#$LumA
#              [,1]     [,2]  [,3]
#Min.    -1.320e+05 -8.23700 19.16
#Max.     1.126e+06  8.39600 43.04
#$Basal
#              [,1]    [,2]  [,3]
#Min.    -23010.000 -8.5370 19.90
#Max.    799300.000  8.2840 40.73
#$LumB
#              [,1]     [,2]  [,3]
#Min.    -4.101e+04 -8.40200 18.78
#Max.     1.022e+06  8.45600 42.52
#$Her2
#              [,1]     [,2]  [,3]
#Min.    -2.515e+04 -8.37300 19.80
#Max.     1.318e+06  8.35300 41.41
#$normal
#              [,1]     [,2]  [,3]
#Min.      -9339.00 -8.42100 20.56
#Max.    1001000.00  8.40100 40.92

#como los subtipos varian en rangos muy diferentes, hay que escalar
concatenadas=lapply(concatenadas,function(x) scale(do.call(rbind,x)))
lapply(concatenadas,function(x) summary(as.numeric(x)))
#$LumA
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-58.2200  -0.0072  -0.0054   0.0000  -0.0043 502.9000 
#$Basal
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-14.9400  -0.0081  -0.0054   0.0000  -0.0037 542.3000 
#$LumB
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-21.0000  -0.0074  -0.0054   0.0000  -0.0040 527.9000 
#$Her2
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-10.7500  -0.0062  -0.0045   0.0000  -0.0035 565.8000 
#$normal
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -4.8530  -0.0071  -0.0050   0.0000  -0.0036 521.0000 
write.table(concatenadas$Her2,"parallel-aracne/data/her2.txt",sep='\t',quote=F)

subtipos=subtipos[subtipos$sample%in%substr(unlist(sapply(concatenadas,colnames)),1,16),]
mir=grep("hsa",rownames(concatenadas$Her2))
mrna=grep("ENSG",rownames(concatenadas$Her2))
concatenadas=lapply(concatenadas,function(x)
 list(mir=x[mir,],mrna=x[mrna,],methy=x[!(1:nrow(x)%in%mir|1:nrow(x)%in%mrna),]))
save(concatenadas,subtipos,file="resultados/porSubtiScaled.RData")

library(parallel)
library(FactoMineR)
library(factoextra)

mir=do.call(cbind,sapply(concatenadas,function(x) x[1]))
mrna=do.call(cbind,sapply(concatenadas,function(x) x[2]))
methy=do.call(cbind,sapply(concatenadas,function(x) x[3]))
omics=list(methy=t(methy),mir=t(mir),mrna=t(mrna))
subtipos=subtipos[order(match(subtipos$sample,substr(rownames(omics$mir),1,16))),]
cl <- makeCluster(7)
clusterEvalQ(cl,{library(FactoMineR)})
PCAomics=parLapply(cl, omics,function(x) PCA(x,ncp=5,graph=F))
stopCluster(cl)
pdf("PCAomics.pdf")
fviz_pca_ind(PCAomics$mrna,geom="point",habillage=subtipos$pbcmc2,palette=ggsci::pal_lancet("lanonc")(5),title="Expresión de mRNAs",pch=19)
fviz_pca_ind(PCAomics$mir,geom="point",habillage=subtipos$pbcmc2,palette=ggsci::pal_lancet("lanonc")(5),title="Expresión de miRNAs",pch=19)
fviz_pca_ind(PCAomics$methy,geom="point",habillage=subtipos$pbcmc2,palette=ggsci::pal_lancet("lanonc")(5),title="Metilación de CpGs",pch=19)
dev.off()
