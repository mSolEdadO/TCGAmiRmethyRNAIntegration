library(data.table)
library(FactoMineR)
library(factorextra)


subtypes=read.table("parallel-aracne/subtype.tsv",header=T,sep='\t')
expre=fread("parallel-aracne/expreNormi.tsv")
miR=fread("parallel-aracne/miRNormi.tsv")
methy=fread("methyM.tsv")
expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)

nrow(expre)
#[1] 17359
nrow(miR)
#[1] 669
nrow(methy)
#[1] 383408
#choose methy order
subtypes=subtypes[order(match(subtypes$barcode,colnames(methy))),]
expre=expre[,order(match(colnames(expre),subtypes$barcode))]
miR=miR[,order(match(colnames(miR),subtypes$barcode))]
names(concatenated)=gsub("BRCA.","",levels(subtypes$subtype))

#data per subtype
concatenated=lapply(levels(subtypes$subtype),function(x) 
	list(methy[,subtypes$subtype==x],expre[,subtypes$subtype==x],
		miR[,subtypes$subtype==x]))

#pcs per subtype for transcripts
PCAomics=lapply(concatenated,function(x) 
	PCA(t(x[[2]]),scale.unit = TRUE,graph=F))
sapply(1:5,function(x) sum(PCAomics[,x]$eig[,3]<75)+1)#PCs needed
#[1]  68  28 169  79  38
sapply(1:5,function(x) length(PCAomics[,x]$eig[,3]))#total PCs
#[1] 125  45 422 145 100
#pcs per subtype for miRNAs
PCAomics=lapply(concatenated,function(x) 
	PCA(t(x[[2]]),scale.unit = TRUE,graph=F))
sapply(PCAomics,function(x) sum(x$eig[,3]<75)+1)
# Basal   Her2   LumA   LumB Normal 
#    64     28    148     74     54 
sapply(PCAomics,function(x) nrow(x$eig))
# Basal   Her2   LumA   LumB Normal 
#   125     45    422    145    100 
pdf("miR_PCA.pdf")
lapply(PCAomics,function(x) fviz_eig(x, addlabels = TRUE))
dev.off()

#matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
sapply(concatenated,dim)
#      Basal   Her2   LumA   LumB Normal
#[1,] 401436 401436 401436 401436 401436
#[2,]    126     46    423    146    101
lapply(1:5,function(x) write.table(concatenated[[x]],
	paste(names(concatenated)[x],"txt",sep='.'),sep='\t',quote=F))


#########################################to go back
files=list.files()
files=files[grep("txt",files)]
library(data.table)
concatenated=lapply(files,fread)
concatenated=lapply(concatenated,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(concatenated)=gsub(".txt","",files)
library(FactoMineR)
concatenated=lapply(concatenated,function(y) apply(cbind(c(1,383409,400768),c(383408,400767,401436)),1,function(x) y[x[1]:x[2],]))
PCAomics=pbapply::pblapply(concatenated,function(x)
 sapply(x,function(y) PCA(t(y),
 scale.unit = TRUE,
 graph=F)))




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
