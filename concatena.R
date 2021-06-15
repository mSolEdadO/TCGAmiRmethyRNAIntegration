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

#to go back
#apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
#	function(x) data[x[1]:x[2],])

#########################################PCs per subtype & data
library(FactoMineR)
library(factoextra)
library(ggplot2)

#function to wrap it all for the different subtypes
check_var=function(data,name){
#data has features in rows & samples in columns
	print(name)
	mfa=MFA(t(data),group=c(393132,17077,604),#size of categories
	#mfa=MFA(t(data),group=c(10,10,10),#size of categories
		name.group=c("CpGs","transcripts","miRNAs"),graph=F,ncp=3)
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
pblapply(1:4,function(x) 
	check_var(concatenated[[x]],names(concatenated)[x]))
#~02d 23h 33m 41s
#Basal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#   13    32    32 
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#58.33845 36.82078 35.84572 
#Weights per omic
#[1] 2.435986e-05 1.842977e-03 6.675847e-02
#Her2
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    15    16
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#Weights per omic
#[1] 2.010012e-05 0.001200744 0.0335368
#LumA
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    20    63    64
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#50.15342 25.43620 21.87975 
#Weights per omic
#[1] 2.799981e-05 1.915935e-03 1.069552e-01
#LumB
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    12    38    36
#
#Normal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    3    12    20
#Variance when 20 PCs are kept
#       CpGs transcripts      miRNAs 
#   78.36429    66.41633    50.37623 
#Weights per omic
#[1] 8.721703e-06 7.684326e-04 4.049163e-02

##################################drop near zero var features???????''

