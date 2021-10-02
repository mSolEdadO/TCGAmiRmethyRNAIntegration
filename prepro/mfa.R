#!/usr/bin/env Rscript
library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)

name = commandArgs(trailingOnly=TRUE)

check_var=function(data,name){
#data has features in rows & samples in columns
	print(name)
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
	write.table(mfa$summary.quanti,paste(name,"ecart",sep='.'),
		sep='\t',quote=F)
	p=ggplot(mfa$summary.quanti,aes(y=ecart.type,
		x=names(mfa$separate.analyses)[mfa$summary.quanti$group]))
	p+geom_boxplot()+ylab("sd")+xlab("")+
		scale_y_continuous(trans="log10")+
	 	theme(text=element_text(size=18))
	ggsave(file=paste(name,"sd.png",sep='.'))
	#mixed PCs
	png(paste(name,"global","png",sep='.'))
	print({
	fviz_screeplot(mfa,addlabels=F,ncp=45,main="global")})
	dev.off()
}

mtrx=fread(paste(name,"mtrx",sep='.'))
mtrx=as.matrix(mtrx[,2:ncol(mtrx)],rownames=mtrx$V1)
check_var(mtrx,name)

##################################RESULTS##############################
#Basal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#   13    32    32 
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#58.33845 36.82078 35.84572 
#Her2
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    8    15    16
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#76.37429 61.49526 61.38255
#LumA
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    20    63    64
#Variance when 20 PCs are kept
#   methy      RNA    miRNA 
#50.15342 25.43620 21.87975 
#LumB
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    12    38    36
#"Variance when 20 PCs are kept"
#       CpGs transcripts      miRNAs 
#   58.67621    32.26742    32.46045 
#Normal
#PCs to keep 50% of variance
#methy   RNA miRNA 
#    3    12    20
#Variance when 20 PCs are kept
#       CpGs transcripts      miRNAs 
#   78.36429    66.41633    50.37623 

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


