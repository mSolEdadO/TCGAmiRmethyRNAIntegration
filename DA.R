library(limma)
#library(sva)
library(missMethyl)
load("porSubti.RData")
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://f1000research.com/articles/5-1408/v1

methy=do.call(cbind,sapply(concatenadas,function(x) x[[2]]))
design=methyDesign[methyDesign$barcode%in%colnames(methy),]
#add gender factor, to control for it on DM, nor on DE coz ARSyN is expected to have wiped all unwanted signal
table(methyDesign[,c(4,10)])
#        gender
#subtype  female male
#  Basal     135    0
#  Her2       74    1
#  LumA      331    0
#  LumB      170    7
#  normal     96    0
load("subtiTMMArsyn.RData")
gender=sapply(as.character(design$patient),function(x) unique(subtipos$gender[as.character(subtipos$patient)==x]))
design=cbind(design,gender)

#set comparisons
DE.design=model.matrix(~0+design$subtype)#a model where each coefficient corresponds to a group mean
colnames(DE.design)=gsub("design.","",colnames(DE.design))
contr.mtrx=makeContrasts(
	basal_normal=subtypeBasal-subtypenormal,
	her2_normal=subtypeHer2-subtypenormal,
	luma_normal=subtypeLumA-subtypenormal,
	lumb_normal=subtypeLumB-subtypenormal,
levels=DE.design)

#mRNA DE
rna=do.call(cbind,sapply(concatenadas,function(x) x[[3]]))
v=voom(rna,DE.design,plot=T,save.plot=T)#ideally count&variance are indi, but a smooth fitted curve is good enough to discard voom 
#counts are more variable at lower expression, voom address this by making the data “normal enough”
#rna variance~[0.1,0.3]
fit=lmFit(rna,DE.design)#fit a linear model using weighted least squares for each gene
fitSubtype = contrasts.fit(fit, contr.mtrx)


#esta bien usar t-test? no sabes si la expresión tiene distro normal, pero seguro no
#tampoco esperas homogeneity no?


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

#miR DE
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

#CpG DM reccomended method: RUV-inverse outperforms existing methods in DA of 450k data
M=lapply(1:4,function(x) cbind(concatenadas[[x]][[2]],concatenadas$normal[[2]]))
i=sapply(M,ncol)-75
DE.design=lapply(1:4,function(x) 
	model.matrix(~1+as.factor(c(rep(names(i)[x],i[x]),rep("normal",75)))))
#stage 1: find empirical controls not associated with subtypes 
fit = lapply(1:4,function(x) lmFit(M[[x]],DE.design[[x]]))
#unadjusted limma is the worst method coz it has lots of false positive [Maksimovic2015]
efit=lapply(fit,eBayes)
topr=lapply(efit,function(x) topTable(x,num=Inf))## Removing intercept from test coefficients
sapply(topr,function(x) sum(x$logFC[x$adj.P.Val<0.01]<0))#most DMs possible
#[1] 118260  82020 113250 109812
sapply(topr,function(x) sum(x$logFC[x$adj.P.Val<0.01]>0))
#[1] 113345 131402 103097 102783
ctl=lapply(topr,function(x) x$adj.P.Val>0.5) 
# performance is consistent when controls are selected based on FDR cut-off
sapply(ctl,table)#Maksimovic uses  2051 ‘true’ positives and 170 629 ‘true’ negatives
#        [,1]   [,2]   [,3]   [,4]
#FALSE 342764 336375 330017 332717
#TRUE   41811  48200  54558  51858
# Stage 2 analysis
rfit <- sapply(1:4,function(x) RUVfit(data=M[[x]], design=DE.design[[x]], coef=2, ctl=ctl[[x]])) 
rfit <- lapply(rfit,RUVadj)
DM.cpg=lapply(rfit,function(x) topRUV(x,num=Inf))
names(DM.cpg1)=c("luma_normal","basal_normal","lumb_normal","her2_normal")
sapply(DM.cpg,function(x) sum(x[x$p.ebayes.BH<0.01,1]<0))
#[1]    8 1146 1710 3620
sapply(DM.cpg,function(x) sum(x[x$p.ebayes.BH<0.01,2]>0))
#[1]   48 1312 1664 2877
#lumA may have less DM coz it has more samples, so the same threshold doesn't work for all subtypes 

#Alternative method
#given the large N(>>5) and cancer effect (gives the 1st PC), all methods are expected to work similarly
#though this only controls for KNOWN batches
DE.design=model.matrix(~0+design$subtype+design$plate+design$sample+design$vial+design$portion,design$gender)
colnames(DE.design)=gsub("design.","",colnames(DE.design))
fit1 = lmFit(methy,DE.design)
contr.mtrx=rbind(contr.mtrx,matrix(rep(0,(ncol(DE.design)-nrow(contr.mtrx))*ncol(contr.mtrx)),ncol=ncol(contr.mtrx)))
#0 for covariates, we're only after subtype signal 
fitSubtype.M1 = contrasts.fit(fit1, contr.mtrx)
tfitSubtype.M1 = treat(fitSubtype.M1,lfc=log2(1.5))
DM.cpg1=lapply(1:4,function(x) topTreat(tfitSubtype.M1,coef=x,n=nrow(methy)))
names(DM.cpg1)=colnames(contr.mtrx)
sapply(DM.cpg1,function(x) sum(x$logFC[x$adj.P.Val<0.01]<0))
#[1] 11266 11289  6882 25298
sapply(DM.cpg1,function(x) sum(x$logFC[x$adj.P.Val<0.01]>0))
#[1] 18367 27116 10847 39200
png("DM.png")
par(mfrow=c(3,4))
sapply(efit,limma::plotMA)
sapply(DM.cpg,function(x) plot(x[,c(1,5)],log="y",ylim=c(1,1e-13),ylab="log(p.ebayes)",xlab="lfc",pch='.'))
sapply(1:4,function(x) limma::volcanoplot(tfitSubtype.M1,coef=x,main=names(DM.cpg1)[x]))
dev.off()

save(DE.genes,DE.miR,DM.cpg,DM.cpg1,file="DA.RData")
