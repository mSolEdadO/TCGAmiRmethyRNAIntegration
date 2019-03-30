library(limma)
load("porSubti.RData")
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://f1000research.com/articles/5-1408/v1

design=as.data.frame(do.call(rbind,sapply(1:5,function(x) cbind(colnames(concatenadas[[x]][[1]]),names(concatenadas)[x]))))
colnames(design)=c("sample","subtype")
DE.design=model.matrix(~0+design$subtype)#a model where each coefficient corresponds to a group mean
colnames(DE.design)=gsub("design.","",colnames(DE.design))
contr.mtrx=makeContrasts(
	basal_normal=subtypeBasal-subtypenormal,
	her2_normal=subtypeHer2-subtypenormal,
	luma_normal=subtypeLumA-subtypenormal,
	lumb_normal=subtypeLumB-subtypenormal,
levels=DE.design)

rna=do.call(cbind,sapply(concatenadas,function(x) x[[3]]))
v=voom(rna,DE.design,plot=T,save.plot=T)#counts are more variable at lower expression, voom address this by making the data “normal enough”
#ideally count&variance are indi, but a smooth fitted curve is good enough to discard voom
#rna variance~[0.1,0.3]
fit=lmFit(rna,DE.design)#fit a linear model using weighted least squares for each gene
fitSubtype = contrasts.fit(fit, contr.mtrx)
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))#fc thresholding + p.val thresholding augments false positives
#lfc=log2(1.2)or log2(1.5)will usually cause most differentially expressed genes to fc => 2-fold,
# depending on the sample size and precision of the experiment
DE.genes=lapply(1:4,function(x) topTreat(tfitSubtype,coef=x,n=nrow(rna)))
sapply(1:4,function(x) sum(DE.genes[[x]]$adj.P.Val<0.01))
#[1] 4904 4357 3452 4722
pdf("DEgenes.pdf")
par(mfrow=c(2,2))
sapply(1:4,function(x) plotMA(fitSubtype,coef=x))
sapply(1:4,function(x) {
	volcanoplot(tfitSubtype,coef=x)
	abline(h=-log2(0.01),col="red")
})
dev.off()

mirna=do.call(cbind,sapply(concatenadas,function(x) x[[1]]))
v.miR=voom(mirna,DE.design,plot=T,save.plot=T)
#mirna variance~[0.3,1.6]
fit=lmFit(v.miR,DE.design)
	fiTemp=lmFit(mirna,DE.design)#pruebame		   
fitSubtype.miR = contrasts.fit(fit, contr.mtrx)
	fiTemp=contrasts.fit(fiTemp, contr.mtrx)#pruebame
tfitSubtype.miR=treat(fitSubtype.miR, lfc = log2(1.1))
DE.miR=lapply(1:4,function(x) topTreat(tfitSubtype.miR,coef=x,n=nrow(mirna)))
sapply(1:4,function(x) sum(DE.miR[[x]]$adj.P.Val<0.01))
#[1] 328 278 259 291
#no voom gives
#[1] 298 259 248 282
pdf("DEmiRs.pdf")
par(mfrow=c(2,2))     
###plotMA antes de voom!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       sapply(1:4,function(x) plotMA(fiTemp,coef=x))#pruebame
sapply(1:4,function(x) plotMA(fitSubtype.miR,coef=x))
sapply(1:4,function(x) {
	volcanoplot(tfitSubtype.miR,coef=x)
	abline(h=-log2(0.01),col="red")
})
dev.off()


	      
#plot var vs mean a mano????????????????????????????????????????????
methy=do.call(cbind,sapply(concatenadas,function(x) x[[2]]))
v.methy=voom(methy,DE.design,plot=T,save.plot=T)
fit=lmFit(methy,DE.design)
fitSubtype.M = contrasts.fit(fit, contr.mtrx)
###
fitSubtype.M = eBayes(fitSubtype.M)#coefficientes has the log fold changes
#DEM=decideTests(fitSubtype.M,method="separate",p.value=3.158632e-08,adjust.method="fdr")
#summary(DEM)
#       basal_normal her2_normal luma_normal lumb_normal
#Down          52119       35019       47546       60853
#NotSig       312330      309121      291504      263110
#Up            31292       51601       56691       71778

svg("voom.svg")
par(mfrow=c(2,1))
plot(v$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="RNA")
lines(v$voom.line,col="red")
plot(v.miR$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="miRNA")
lines(v.miR$voom.line,col="red")#se espera un curva suave y concava y eso es convexo
#esto ruega que limpie bien los mirna????
dev.off()


#saco los elementos más diferentes (y signif) entre tej normal y cancer
fitSubtype$fdr=apply(fitSubtype$"p.value",2,p.adjust,"fdr")
g=names(which(rowSums(abs(fitSubtype$coefficients)>2)>0
	&rowSums(abs(fitSubtype$fdr)<0.05)>0))
#F-statistic tests whether any of the contrasts are non-zero. With many contrasts, it may be desirable to select genes firstly on the basis of their moderated F-statistics
fitSubtype.miR$fdr=apply(fitSubtype.miR$"p.value",2,p.adjust,"fdr")
mi=names(which(rowSums(abs(fitSubtype.miR$coefficients)>2)>0
	&rowSums(abs(fitSubtype.miR$fdr)<0.05)>0))
fitSubtype.M$fdr=apply(fitSubtype.M$"p.value",2,p.adjust,"fdr")
M=names(which(rowSums(abs(fitSubtype.M$coefficients)>2)>0
	&rowSums(abs(fitSubtype.M$fdr)<0.05)>0))
save(fitSubtype,fitSubtype.miR,fitSubtype.M,file="diff/contrsFiteBay.RData")

rna=rna[rownames(rna)%in%g,]
mirna=mirna[rownames(mirna)%in%mi,]
methy=methy[rownames(methy)%in%M,]

#save(rna,mirna,methy,design,file="conca/subsetmasdiff.Rda")#lo borre porque ya tengo los subsets en useR/MFA.Rda:subti
