library(data.table)
library(FactoMineR)
library(factorextra)


subtype=read.table("subtype.tsv",header=T,sep='\t')
expre=fread("parallel-aracne/RNAseqnormalized.tsv")
miR=fread("parallel-aracne/miRNAseqNormi.tsv")
methy=fread("methyM.tsv")
expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
#choose methy order
subtype=subtype[order(match(subtype$samples,colnames(methy))),]
expre=expre[,order(match(colnames(expre),subtype$samples))]
miR=miR[,order(match(colnames(miR),subtype$samples))]

#data per subtype
concatenated=lapply(levels(subtype$subtype),function(x) 
	list(methy=methy[,subtype$subtype==x],
		RNA=expre[,subtype$subtype==x],
		miRNA=miR[,subtype$subtype==x]))
names(concatenated)=levels(subtype$subtype)

#########################################PCs per subtype & data
her2MFA=MFA(t(her2),group=c(393132,17077,604),#size of categories
	name.group=c("methy","RNA","miRNA"),graph=F,ncp=3)
#PCs to keep 50% of variance
sapply(her2MFA$separate.analyses,function(x) sum(x$eig[,3]<50))+1
#methy   RNA miRNA 
#    8    15    16
pdf("PCsHer2.pdf")
lapply(1:3,function(x) fviz_eig(her2MFA$separate.analyses[[x]],
	addlabels=F,ncp=45,main=names(her2MFA$separate.analyses)[x]))
dev.off()
#what if I keep 20 PC
sapply(her2MFA$separate.analyses,function(x) x$eig[20,3])#var explained
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255 
#mixed PCs
png("PCsHer2.png")
 fviz_screeplot(her2MFA,addlabels=F,ncp=45,main="global")
dev.off()


#For transcripts
PCArna=lapply(concatenated,function(x) 
	PCA(t(x[[2]]),scale.unit = TRUE,graph=F))
sapply(PCArna,function(x) sum(x$eig[,3]<50)+1)
#    32     15     63     38     12 
pdf("PCrna.pdf")
lapply(1:5,function(x) fviz_eig(PCArna[[x]],addlabels=F,
	ncp=50,main=names(PCArna)[x]))
dev.off()
#what if I keep 25 PC
sapply(PCArna,function(x) x$eig[25,3])#var explained
#   Basal     Her2     LumA     LumB   Normal 
#42.90686 71.39674 29.40306 37.90650 72.93814 

#For miRNAs
PCAmir=lapply(concatenated,function(x) 
	PCA(t(x[[3]]),scale.unit = TRUE,graph=F))
sapply(PCAmir,function(x) sum(x$eig[,3]<50)+1)
# Basal   Her2   LumA   LumB Normal 
#    32     16     64     36     20 
pdf("PCmiRNA.pdf")
lapply(1:5,function(x) fviz_eig(PCAmir[[x]],addlabels=F,
	ncp=50,main=names(PCAmir)[x]))
dev.off()
#what if I keep 25 PC
sapply(PCAmir,function(x) x$eig[25,3])#var explained
#   Basal     Her2     LumA     LumB   Normal 
#42.28663 71.59183 26.04078 38.46184 58.49842 
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
 cbind(names(PCAmir)[x],
 PCAmir[[x]]$call$ecart.type))))
colnames(temp)=c("subtype","sde")
temp$sde=as.numeric(as.character(temp$sde))
png("standardErrormiR.png")
ggplot(temp,aes(x=sde))+
geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)
dev.off()










#########################################drop near zero var features

##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
sapply(concatenated,dim)
#      Basal   Her2   LumA   LumB Normal
#[1,] 410813 410813 410813 410813 410813
#[2,]    128     46    416    140     75
lapply(1:5,function(x) write.table(concatenated[[x]],
	paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))

#########################################to go back
files=list.files()
files=files[grep("txt",files)]
library(data.table)
concatenated=lapply(files,fread)
concatenated=lapply(concatenated,function(x) 
	as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(concatenated)=gsub(".txt","",files)
concatenated=lapply(concatenated,function(y) 
	apply(cbind(c(1,383409,400768),c(383408,400767,401436)),1,
		function(x) y[x[1]:x[2],]))
library(FactoMineR)
PCAomics=pbapply::pblapply(concatenated,function(x)
 sapply(x,function(y) PCA(t(y),
 scale.unit = TRUE,
 graph=F)))

#################################
###################3

#lapply(concatenadas,function(x) sapply(x,function(y) summary(as.numeric(y))))
#$normal
#           methy transcri        mir
#Min.    -8.42100 -39100.0  -7370.000
#Max.     8.40100 370900.0 911600.000
#$LumA
#           methy   transcri        mir
#Min.    -8.37300 -1.046e+06 -1.646e+05
#Max.     8.45200  3.383e+06  9.950e+05
#$LumB
#           methy   transcri        mir
#Min.    -8.40200 -1.130e+06 -19390.000
#Max.     8.45600  2.446e+06 937100.000
#$Her2
#          methy   transcri        mir
#Min.    -8.3730 -719500.00 -1.864e+04
#Max.     8.3530 3143000.00  1.068e+06
#$Basal
#          methy   transcri        mir
#Min.    -8.5370 -1.196e+06 -26820.000
#Max.     8.2840  3.294e+06 680800.000

#since omics have different ranges →MFA weight
evals=lapply(concatenadas,function(x) sapply(x,function(y) svd(y,nu=0,nv=0)))
evals=lapply(evals,function(x) sapply(x,function(y) 1/sqrt(y[1])))
#$LumA
#     methy.d   transcri.d        mir.d 
#0.0049129995 0.0002105128 0.0002012404 
#$LumB
#     methy.d   transcri.d        mir.d 
#0.0066157577 0.0002897397 0.0002886322 
#$normal
#     methy.d   transcri.d        mir.d 
#0.0072579414 0.0003289092 0.0003312752 
#$Her2
#     methy.d   transcri.d        mir.d 
#0.0085562074 0.0003376987 0.0003486907 
#$Basal
#     methy.d   transcri.d        mir.d 
#0.0065336421 0.0002825550 0.0003301057 
concatenadas=lapply(1:5,function(x) scale(do.call(rbind,lapply(1:3,function(y) 
	as.matrix(concatenadas[[x]][[y]]*evals[[x]][y])))))
save(concatenadas,file="escaladas.RData")
lapply(concatenadas,function(x) summary(as.numeric(x)))
#$LumA
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-49.0892  -0.0787  -0.0445   0.0000  -0.0225 396.1100 
#$LumB
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -6.3675  -0.0844  -0.0489   0.0000  -0.0257 396.6516 
#normal
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2.9152  -0.0843  -0.0471   0.0000  -0.0237 379.8674 
#Her2
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -6.0729  -0.0750  -0.0427   0.0000  -0.0227 413.7628 
#Basal
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-13.1357  -0.0875  -0.0497   0.0000  -0.0257 369.4758 
