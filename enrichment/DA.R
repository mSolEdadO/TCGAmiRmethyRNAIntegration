library(limma)
library(data.table)
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://f1000research.com/articles/5-1408/v1

#################RNA###########################################
subtype=read.table("../Downloads/subtype.tsv",header=T,sep='\t')
expre=fread("../Downloads/RNAseqnormalized.tsv")
expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)

#set comparisons
#~0 gives a model where each coefficient corresponds to a group mean
subtype$subtype=factor(subtype$subtype)
design=model.matrix(~0+subtype,subtype)
#fix names
colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
	Basal_Normal=Basal-Normal,
	Her2_Normal=Her2-Normal,
	LumA_Normal=LumA-Normal,
	LumB_Normal=LumB-Normal,
levels=design)

#check if count&variance are indi
#if counts are more variable at lower expression, voom makes the 
#data “normal enough”, 
#wont work with ur log2-normalized stuff
#v=voom(expre,design,plot=T,save.plot=T)#no need if fitted curve is smooth 
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

#fit a linear model using weighted least squares for each gene
fit=lmFit(log2(expre),design)#lots of NA, if no log2 lfc is huge & u get lot of DE.genes
fitSubtype = contrasts.fit(fit, contr.mtrx)
#treat is better than fc+p.val thresholds, that increase FP
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))
#log2(1.2 or 1.5) will usually give DE genes with fc => 2
#depending on the sample size and precision of the experiment
DE.genes=lapply(1:4,function(x) 
	topTreat(tfitSubtype,coef=x,n=nrow(expre)))
names(DE.genes)=colnames(contr.mtrx)
sapply(DE.genes,function(x) sum(x$adj.P.Val<0.01))
#        5334         4455         4075         4982 
##pdf("DEgenes.pdf")
#par(mfrow=c(2,2))
#sapply(1:4,function(x) plotMA(fitSubtype,coef=x))
#sapply(1:4,function(x) {
#	volcanoplot(tfitSubtype,coef=x)
#	abline(h=-log2(1.2),col="red")
#})
#dev.off()
temp=do.call(rbind,lapply(1:4,function(x) 
	cbind(contrast=names(DE.genes)[x],
		ensembl_gene_id=rownames(DE.genes[[x]]),
		DE.genes[[x]])))
png("logFC.png")
ggplot(temp,aes(y=logFC,x=contrast,color=contrast))+
	geom_boxplot()+theme(legend.position="none")
dev.off()
write.table(temp,"DE.genes.tsv",sep='\t',quote=F,row.names=F)
#next:GSEA
#################miRNA###########################################
mir=fread("../Downloads/miRNAseqNormi.tsv")
mir=as.matrix(mir[,2:ncol(mir)],rownames=mir$V1)

v=voom(mir,design,plot=T,save.plot=T)#coz mir is normalized, but no log2 transformed
fit=lmFit(v,design)
fitSubtype = contrasts.fit(fit, contr.mtrx)
#treat doesn't find any DE miR
#tfitSubtype=treat(fitSubtype, lfc = log2(1.2))
#DE.miR=lapply(1:4,function(x) 
#	topTreat(tfitSubtype,coef=x,n=nrow(mir)))
temp=eBayes(fitSubtype)
DE.miR=lapply(1:4,function(x) 
	topTable(temp,coef=x,n=Inf))
names(DE.miR)=colnames(contr.mtrx)
sapply(DE.miR,function(x) sum(x$adj.P.Val<0.05))
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#         174            0           39           75 
temp=do.call(rbind,lapply(1:4,function(x) 
	cbind(contrast=names(DE.miR)[x],
		id=rownames(DE.miR[[x]]),
		DE.miR[[x]])))
write.table(temp,"DE.miR.tsv",sep='\t',quote=F,row.names=F)

#################methylation###########################################
clin <- TCGAbiolinks::GDCquery_clinic("TCGA-BRCA","clinical")
colnames(subtype)[3]="bcr_patient_barcode"
subtype=merge(subtype,clin,by="bcr_patient_barcode",all.x=T)
#apply(subtype,2,table) choose cols with various categories
subtype=subtype[,c("bcr_patient_barcode",
					"samples",
					"tissue",
					"subtype",
					"age_at_diagnosis",
					"treatments_pharmaceutical_treatment_or_therapy",
					"treatments_radiation_treatment_or_therapy",
					"race",
					"ajcc_pathologic_m",
					"ajcc_pathologic_n",
					"morphology",
					"ajcc_pathologic_t",
					"prior_malignancy",
					"primary_diagnosis",
					"ajcc_pathologic_stage")]
write.table(subtype,"../Downloads/subtype.tsv",sep='\t',quote=F,row.names=F)

methy=data.table::fread("methyM.tsv")
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#change category names for numbers or makeContrast wont work 
subtype=cbind(subtype[,1:5],
			  apply(subtype[,6:15],2,function(x)
 			     factor(as.numeric(factor(x)))))
subtype[is.na(subtype)]=0#force NA to be another category or ur loose samples

#given the large N(>>5) and cancer effect, all methods are expected to work similarly
design=model.matrix(~0+subtype+age_at_diagnosis+morphology+
	ajcc_pathologic_m+ajcc_pathologic_n+ajcc_pathologic_t+
	primary_diagnosis+ajcc_pathologic_stage,subtype)
#used this complex model coz methy didn't went through ARSYN 
#no treatment covarites/prior_malignancy coz groups're mixed in PCA
colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
	Basal_Normal=Basal-Normal,
	Her2_Normal=Her2-Normal,
	LumA_Normal=LumA-Normal,
	LumB_Normal=LumB-Normal,
levels=design)

fit=lmFit(methy,design)
#Coefficients not estimable: primary_diagnosis
#Partial NA coefficients for 393132 probe(s) 
fitSubtype = contrasts.fit(fit, contr.mtrx)
#temp=eBayes(fitSubtype) #only with subtype
#DMcpgs=lapply(1:4,function(x) topTable(temp,coef=x,n=Inf))
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#      246810       228203       261208       256139 
#use treat to get lower dm cpgs
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))
DMcpgs=lapply(1:4,function(x) 
	topTreat(tfitSubtype,coef=x,n=Inf))
names(DMcpgs)=colnames(contr.mtrx)
sapply(DMcpgs,function(x) sum(x$adj.P.Val<0.05))
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#       61132        77776        89863       114345 
#when desing only considers subtype, DMs don't change much
#Maksimovic2015 says unadjusted limma has lots of FP
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#       64687        79396        88830       115412 
temp=do.call(rbind,lapply(1:4,function(x) 
	cbind(contrast=names(DMcpgs)[x],
		id=rownames(DMcpgs[[x]]),
		DMcpgs[[x]])))
write.table(temp,"DMcpgs.tsv",sep='\t',quote=F,row.names=F)



#RUV-inverse outperforms existing methods in DA of 450k data
design=model.matrix(~0+subtype,subtype)
colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
	Basal_Normal=Basal-Normal,
	Her2_Normal=Her2-Normal,
	LumA_Normal=LumA-Normal,
	LumB_Normal=LumB-Normal,
levels=design)
#stage 1: find empirical controls not associated with subtypes 
fit=lmFit(methy,design)
fitSubtype = contrasts.fit(fit, contr.mtrx)
temp=eBayes(fitSubtype) #only with subtype
DMcpgs=lapply(1:4,function(x) topTable(temp,coef=x,n=Inf))
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#      246810       228203       261208       256139 
ctlimma=lapply(DMcpgs,function(x) x$adj.P.Val > 0.5)
sapply(ctlimma,table)
#      Basal_Normal Her2_Normal LumA_Normal LumB_Normal
#FALSE       339154      329948      345252      338195
#TRUE         53978       63184       47880       54937

# Stage 2 analysis
library(missMethyl)
rfit <- sapply(1:4,function(x) 
	RUVfit(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(x,5)]],
		   X=factor(subtype$subtype[subtype$subtype%in%unique(subtype$subtype)[c(x,5)]]),
		   ctl=ctlimma[[x]])) 
rfit <- lapply(rfit,RUVadj)
DMcpg=lapply(rfit,function(x) topRUV(x,num=Inf))
names(DMcpgs)=colnames(contr.mtrx)
sapply(DMcpg,function(x) sum(x[x$p.ebayes.BH<0.01,1]<0))

methy=fread("methyM.tsv",select=which(subtype$subtype%in%c("Basal","Normal")))



#library(sva)
load("porSubti.RData")

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

#CpG DM reccomended method: RUV-inverse outperforms existing methods in DA of 450k data
M=lapply(1:4,function(x) cbind(concatenadas[[x]][[2]],concatenadas$normal[[2]]))
i=sapply(M,ncol)-75
DE.design=lapply(1:4,function(x) 
	model.matrix(~1+as.factor(c(rep(names(i)[x],i[x]),rep("normal",75)))))
#stage 1: find empirical controls not associated with subtypes 
fit = lapply(1:4,function(x) lmFit(M[[x]],DE.design[[x]]))
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
