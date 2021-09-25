#!/usr/bin/env Rscript
library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)

files=list.files()
files=files[grep("mtrx",files)]
data=lapply(files,fread)
data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
data=do.call(cbind,data)#all subtypes together

mfa=MFA(t(data),group=c(393132,17077,604),#size of categories
	#mfa=MFA(t(data),group=c(10,10,10),#size of categories
		name.group=c("CpGs","transcripts","miRNAs"),graph=F,ncp=3)

print("PCs to keep 50% of variance")
	print(sapply(mfa$separate.analyses,function(x)
	 	sum(x$eig[,3]<50))+1)
print("PCs to keep 75% of variance")
	print(sapply(mfa$separate.analyses,function(x)
	 	sum(x$eig[,3]<75))+1)
pdf(paste(name,"pdf",sep='.'))
print({
lapply(1:3,function(x) fviz_eig(mfa$separate.analyses[[x]],
	addlabels=F,ncp=50,main=names(mfa$separate.analyses)[x]))
})
dev.off()
write.table(mfa$summary.quanti,paste(name,"ecart",sep='.'),sep='\t',quote=F)
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
normalized=data/mfa$global.pca$call$col.w#var=1 for all groups
subtype=read.table("subtype.tsv",header=T,sep='\t')
lapply(levels(subtype$subtype),function(x)
	write.table(normalized[,colnames(normalized)%in%subtype$samples[subtype$subtype==x]],
	paste(x,"normalized",sep='.'),sep='\t',quote=F))
#but sgcca output is the same for scaled data normalized data

##################################RESULTS##############################
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
#"Variance when 20 PCs are kept"
#       CpGs transcripts      miRNAs 
#   58.67621    32.26742    32.46045 
#[1] "Weights per omic"
#[1] 1.547911e-05 2.198907e-03 7.503839e-02
#Normal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    3    12    20
#Variance when 20 PCs are kept
#       CpGs transcripts      miRNAs 
#   78.36429    66.41633    50.37623 
#Weights per omic
#[1] 8.721703e-06 7.684326e-04 4.049163e-02

###############################IN PARALLEL WITH CONDOR
#bash script
#Rscript mfa.R $1

#sub file
# Opciones generales de HTCondor.
#universe = vanilla
#initialdir = /home/msoledad/parallel-R
#should_transfer_files = NO
#getenv = True
# Recursos necesarios para ejecutar el trabajo.
#request_cpus = 1
#request_memory = 3GB
#request_disk = 1GB
####repeat next lines per subtype
#Executable = mfa.sh
#arguments= "LumB"
#Error = mfa.$(Process).err
#Output = mfa.$(Process).out
#Log = mfa.$(Process).log
#Queue


