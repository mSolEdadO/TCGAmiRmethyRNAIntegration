library(limma)
library(sva)
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
v=voom(rna,DE.design,plot=T,save.plot=T)#ideally count&variance are indi, but a smooth fitted curve is good enough to discard voom 
#counts are more variable at lower expression, voom address this by making the data “normal enough”
#rna variance~[0.1,0.3]
fit=lmFit(rna,DE.design)#fit a linear model using weighted least squares for each gene
fitSubtype = contrasts.fit(fit, contr.mtrx)
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))#fc+p.val thresholding increases false positives, treat overpasses this
#lfc=log2(1.2)or log2(1.5)will usually cause most differentially expressed genes to fc => 2-fold,
# depending on the sample size and precision of the experiment
DE.genes=lapply(1:4,function(x) topTreat(tfitSubtype,coef=x,n=nrow(rna)))
names(DE.genes)=c("basal_normal","her2_normal","luma_normal","lumb_normal")
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
v.miR=voom(mirna,DE.design,plot=T,save.plot=T)#mirna variance~[0.3,1.6]
fit=lmFit(v.miR,DE.design)
fitSubtype.miR = contrasts.fit(fit, contr.mtrx)
tfitSubtype.miR=treat(fitSubtype.miR, lfc = log2(1.1))
DE.miR=lapply(1:4,function(x) topTreat(tfitSubtype.miR,coef=x,n=nrow(mirna)))
names(DE.miR)=c("basal_normal","her2_normal","luma_normal","lumb_normal")
sapply(1:4,function(x) sum(DE.miR[[x]]$adj.P.Val<0.01))
#[1] 328 278 259 291
#no voom gives
#[1] 298 259 248 282
fiTemp=lmFit(mirna,DE.design) 
fiTemp=contrasts.fit(fiTemp, contr.mtrx)
pdf("DEmiRs.pdf")
par(mfrow=c(2,2))
sapply(1:4,function(x) plotMA(fiTemp,coef=x))
sapply(1:4,function(x) plotMA(fitSubtype.miR,coef=x))
sapply(1:4,function(x) {
	volcanoplot(tfitSubtype.miR,coef=x)
	abline(h=-log2(0.01),col="red")
})
dev.off()
save(DE.genes,DE.miR,file="DA.RData")

svg("voom.svg")
par(mfrow=c(2,1))
plot(v$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="RNA")
lines(v$voom.line,col="red")
plot(v.miR$voom.xy,xlab="log2(count size + 0.5)",ylab="Sqrt(standard deviation)",main="miRNA")
lines(v.miR$voom.line,col="red")
dev.off()





methy=do.call(cbind,sapply(concatenadas,function(x) x[[2]]))
mod0=model.matrix(~0,data=as.data.frame(t(methy)))
n.sv=num.sv(methy,DE.design)#100
svobj = sva(methy,DE.design,mod0,n.sv=n.sv,B=100)#estimate surrogate variables = estimable 
#linear combinations of the true unmeasured or unmodeled factors causing noise.
modSv = cbind(DE.design,svobj$sv)
fit = lmFit(methy,modSv)
contr.mtrx=rbind(contr.mtrx[1:5,],matrix(rep(0,n.sv*ncol(contr.mtrx)),ncol=4))#0 for sv, just added for coherence with fit than nows has 105 rows
fitSubtype.M = contrasts.fit(fit, contr.mtrx)
tfitSubtype.M = treat(fitSubtype.M,lfc=log2(1.5))
fiTemp=lmFit(methy,DE.design) 
fiTemp=contrasts.fit(fiTemp, contr.mtrx)
fitTemp = treat(fiTemp,lfc=log2(1.5))
 DM.cpg0=lapply(1:4,function(x) topTreat(fitTemp,coef=x,n=nrow(methy)))
DM.cpg=lapply(1:4,function(x) topTreat(tfitSubtype.M,coef=x,n=nrow(methy)))
lapply(DM.cpg,function(x) sum(x$adj.P.Val<0.01))
[[1]]
[1] 1
lapply(DM.cpg0,function(x) sum(x$adj.P.Val<0.01))
[1] 53355

[[2]]
[1] 73661

[[3]]
[1] 72595

[[4]]
[1] 105496

design=cbind(design,sapply(strsplit(as.character(design$sample),"-"),function(x) x[6]))
colnames(design)[3]="plate"
design=cbind(design,sapply(strsplit(as.character(design$sample),"-"),function(x) substr(x[4],1,2)))
colnames(design)[4]="sample"
design=cbind(design,sapply(strsplit(as.character(design$sample),"-"),function(x) substr(x[4],3,3)))
colnames(design)[5]="vial"
design=cbind(design,sapply(strsplit(as.character(design$sample),"-"),function(x) substr(x[5],1,2)))
colnames(design)[6]="portion"
DE.design=model.matrix(~0+design$subtype+design$plate+design$sample+design$vial+design$portion)
colnames(DE.design)=gsub("design.","",colnames(DE.design))
fit1 = lmFit(methy,DE.design)
contr.mtrx=rbind(contr.mtrx[1:5,],matrix(rep(0,61*ncol(contr.mtrx)),ncol=4))#0 for sv, just added for coherence with fit than nows has 105 rows
fitSubtype.M1 = contrasts.fit(fit1, contr.mtrx)
tfitSubtype.M1 = treat(fitSubtype.M1,lfc=log2(1.5))
DM.cpg1=lapply(1:4,function(x) topTreat(tfitSubtype.M1,coef=x,n=nrow(methy)))
lapply(DM.cpg1,function(x) sum(x$adj.P.Val<0.01))
[[1]]
[1] 77526
[[2]]
[1] 99114
[[3]]
[1] 84615
[[4]]
[1] 124304
pdf("DM.png")
par(mfrow=c(3,4))
sapply(1:4,function(x) limma::plotMA(fiTemp,coef=x))
sapply(1:4,function(x) limma::plotMA(fitSubtype.M,coef=x))
sapply(1:4,function(x) limma::plotMA(fitSubtype.M1,coef=x))
dev.off()
#save(rna,mirna,methy,design,file="conca/subsetmasdiff.Rda")#lo borre porque ya tengo los subsets en useR/MFA.Rda:subti
