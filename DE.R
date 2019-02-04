library(limma)
load("conca/porSubti.RData")


rna=do.call(cbind,sapply(concatenadas,function(x) x[[3]]))
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
#tempRNA=voom(rna,DE.design,save.plot=T)#la curva debe ser suave, sino, hay que filtrar más 
fit=lmFit(rna,DE.design)#fit a linear model using weighted least squares for each gene
fitSubtype = contrasts.fit(fit, contr.mtrx)
fitSubtype = eBayes(fitSubtype)
#DEG=decideTests(fitSubtype,method="separate",p.value=1e-100,adjust.method="fdr")#alpha de https://github.com/CSB-IG/tcgarnaseqbc/blob/master/DifGenes.R
#Lower alpha levels are sometimes used when you are carrying out multiple tests at the same time. A common approach is to divide the alpha level by the number of tests being carried out.
#summary(DEG)
#       basal_normal her2_normal luma_normal lumb_normal
#Down            652         426         518         797
#NotSig        12494       13057       12988       12462
#Up              758         421         398         645
#vennDiagram(DEGlimma)

mirna=do.call(cbind,sapply(concatenadas,function(x) x[[1]]))
#tempmiR=voom(mirna,DE.design,save.plot=T)
fit=lmFit(mirna,DE.design)
fitSubtype.miR = contrasts.fit(fit, contr.mtrx)
fitSubtype.miR = eBayes(fitSubtype.miR)
#DEmiR=decideTests(fitSubtype.miR,method="separate",p.value=0.01,adjust.method="fdr")
#alpha de Drago-Garcia
#summary(DEmiR)
#       basal_normal her2_normal luma_normal lumb_normal
#Down             96          44          79          80
#NotSig         1353        1503        1489        1469
#Up              139          41          20          39
methy=do.call(cbind,sapply(concatenadas,function(x) x[[3]]))
fit=lmFit(methy,DE.design)
fitSubtype.M = contrasts.fit(fit, contr.mtrx)
fitSubtype.M = eBayes(fitSubtype.M)
#DEM=decideTests(fitSubtype.M,method="separate",p.value=3.158632e-08,adjust.method="fdr")
#summary(DEM)
#       basal_normal her2_normal luma_normal lumb_normal
#Down          52119       35019       47546       60853
#NotSig       312330      309121      291504      263110
#Up            31292       51601       56691       71778

svg("voom.svg")
par(mfrow=c(2,1))
plot(tempRNA$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="RNA")
lines(tempRNA$voom.line,col="red")
plot(tempmiR$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="miRNA")
lines(tempmiR$voom.line,col="red")#se espera un curva suave y concava y eso es convexo
#esto ruega que limpie bien los mirna????
dev.off()
save(fitSubtype,fitSubtype.miR,fitSubtype.M,file="diff/contrsFiteBay.RData")


#saco los 100 elementos más diferentes entre tej normal y cancer
g=rownames(topTable(fitSubtype,number=100,adjust="fdr",sort.by="F"))#mayor adj.p.val 1.577154e-07
#F-statistic tests whether any of the contrasts are non-zero. With many contrasts, it may be desirable to select genes firstly on the basis of their moderated F-statistics
mi=rownames(topTable(fitSubtype.miR,number=100,adjust="fdr",sort.by="F"))#mayor adj.p.val 2.615420e-186
M=rownames(topTable(fitSubtype.M,number=100,adjust="fdr",sort.by="F"))#mayor adj.p.val 0 
rna=rna[rownames(rna)%in%g,]
mirna=mirna[rownames(mirna)%in%mi,]
methy=methy[rownames(methy)%in%M,]

rna=lapply(unique(design$subtype),function(x) rna[,design$subtype==x])
names(rna)=unique(design$subtype)
mirna=lapply(unique(design$subtype),function(x) mirna[,design$subtype==x])
names(mirna)=unique(design$subtype)
methy=lapply(unique(design$subtype),function(x) methy[,design$subtype==x])
names(methy)=unique(design$subtype)
save(rna,mirna,methy,file="conca/subsetmasdiff.Rda")
