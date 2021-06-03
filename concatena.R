library(data.table)

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

##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
sapply(concatenated,dim)
#      Basal   Her2   LumA   LumB Normal
#[1,] 410813 410813 410813 410813 410813
#[2,]    128     46    416    140     75
lapply(1:5,function(x) write.table(concatenated[[x]],
	paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))

#########################################PCs per subtype & data
library(FactoMineR)
library(factorextra)
library(ggplot2)

#function to wrap it all for the different subtypes
check_var=function(data,name){
#data has features in rows & samples in columns
	print(name)
	mfa=MFA(t(data),group=c(393132,17077,604),#size of categories
	#mfa=MFA(t(data),group=c(10,10,10),#size of categories
		name.group=c("methy","RNA","miRNA"),graph=F,ncp=3)
	print("PCs to keep 50% of variance")
	print(sapply(mfa$separate.analyses,function(x)
	 	sum(x$eig[,3]<50))+1)
	print("Variance when 20 PCs are kept")
	print(sapply(mfa$separate.analyses,function(x) 
		x$eig[20,3]))#var explained
	#	x$eig[5,3]))#var explained
	#elbow plots
	pdf(paste(name,"pdf",sep='.'))
	print({
	lapply(1:3,function(x) fviz_eig(mfa$separate.analyses[[x]],
		addlabels=F,ncp=50,main=names(mfa$separate.analyses)[x]))
	})
	dev.off()
	#sd per omic
	p=ggplot(mfa$summary.quanti,aes(y=ecart.type,
		x=names(mfa$separate.analyses)[mfa$summary.quanti$group]))
	p+geom_boxplot()+ylab("sd")+xlab("")+
		scale_y_continuous(trans="log10")+
	 	theme(text=element_text(size=18))
	ggsave(file=paste(name,"sd.png",sep='.'))
	print("Weights per omic")
	print(unique(mfa$global.pca$call$col.w))
	#mixed PCs
	png(paste(name,"global","png",sep='.'))
	print({
	fviz_screeplot(mfa,addlabels=F,ncp=45,main="global")})
	dev.off()
#var distributions should be similar if the planned 
#approach relies on correlation networks [Tarazona2020]
#normalized=her2/her2MFA$global.pca$call$col.w#var=1 for all groups
#write.table(normalized,"Her2.normalized",sep='\t',quote=F)
#but sgcca output is the same for scaled data normalized data
}

#load a matrix per subtype
files=list.files()
files=files[grep("mtrx",files)]
concatenated=lapply(files,fread)
concatenated=lapply(concatenated,function(x) 
	as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(concatenated)=gsub(".mtrx","",files)
pblapply(1:4,function(x) check_var(concatenated[[x]],names(concatenated)[x]))
#Basal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    32    32
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic

#Her2
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    15    16
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic
# 
#LumA
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    63    64
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic

#LumB
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    38    36
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic

#Normal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    12    20
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic

#########################################drop near zero var features

#########################################to go back
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



##############
#hay que correr lo que sigue ompimizando penalties
 cca=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=10,scale=T)
#da lo mismo normalizada por eigenvalue que no normalizada
#en features seleccionadas y en loadging values
 cca=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=10,scale=F)
#selecciona todos los features de cada bloque

#què significan los links entre nodos????
#tiene caso normalizar???