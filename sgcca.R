#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
subtype=args[1]
BP=args[2]
comps=args[3]

######################MATRIX TO LIST OF MOLECULAR LEVELS
library(data.table)
#data=fread(paste(subtype,"normalized",sep='.'))
data=fread(paste(subtype,"normalized",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
tfs=read.table("TFtargets.tsv",header=T,sep='\t')#https://github.com/slowkow/tftargets
tfs=unique(unlist(strsplit(as.character(tfs$TF),',')))
#length(tfs)
#[1] 991 lot less than https://doi.org/10.1016/j.cell.2018.01.029
#BUT tftargets does have an evidence based list of interactions
data$TFs=data$transcripts[,colnames(data$transcripts)%in%tfs]

#######KEEP ONLY THE TRANSCRIPTS FOR A BIOLOGICAL PROCESS
library(HTSanalyzeR)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    
GS_GO_BP=unlist(GS_GO_BP[names(GS_GO_BP)==BP])
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
	mart=mart)
myannot=myannot[myannot$entrezgene_id%in%unlist(GS_GO_BP),]
data$BP=data$transcripts[,colnames(data$transcripts)%in%
	myannot$ensembl_gene_id[myannot$entrezgene_id%in%GS_GO_BP]])
data$transcripts=NULL

######################CHOOSE SPARSITY VALUES
#sparsity parameters are chosen for each of the 10 MCCV iterations
#using an internal 5-fold CV loop: the parameters that minimize the
#prediction error [Tenenhaus2014]

library(ggplot2)
library(gridExtra)
library(mixOmics)
library(tidyverse)
library(pbapply)

#take model descriptors<-----------------recycled
describe=function(resus,names,omic){
	AVE=as.data.frame(sapply(resus,function(x) 
		do.call(rbind,x$AVE[[1]])))
	nfeat=as.data.frame(sapply(resus,function(x) 
		sapply(x$loadings,function(y) sum(y!=0))))
	#format to plot
	colnames(AVE)=names
	AVE$omic=rownames(nfeat)
	AVE=AVE%>%pivot_longer(-omic,names_to="sparsity",values_to="AVE")
	colnames(nfeat)=names
	nfeat$omic=rownames(nfeat)
	nfeat=nfeat%>%pivot_longer(-omic,names_to="sparsity",
		values_to="nfeatures")
	toplot=merge(nfeat,AVE,by=c("omic","sparsity"))
	#force order
	#toplot$omic=factor(toplot$omic,levels=c("CpGs","transcripts","miRNAs"))
	toplot$sparsity=as.numeric(toplot$sparsity)
	toplot=toplot[toplot$omic==omic,]
return(toplot)}

grid=c(seq(0,.1,.02),seq(0.2,.9,.1))
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(x,1,1,1),scale=T,scheme="centroid"))
#centroid scheme enable two components to be negatively correlated
results=list()
results$CpGs=describe(temp,grid,"CpGs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,x,1,1),scale=T,scheme="centroid"))
results$miRNAs=describe(temp,grid,"miRNAs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,1,x,1),scale=T,scheme="centroid"))
results$TFs=describe(temp,grid,"TFs")
results=do.call(rbind,results)

plots=lapply(levels(results$omic),function(i) 
	ggplot(results[results$omic==i,],aes(x=nfeatures,
		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
		scale_x_continuous(trans="log10")+
		theme(text=element_text(size=18)))
png(paste(subtype,BP,"png",sep='.'))
	print({grid.arrange(plots[[1]],plots[[2]],plots[[3]])})
 #no facet_wrap coz u want indy axes
dev.off()

#turn it analyticalish by choosing the largest fall in AVE
sc=which.max(sapply(14:2,function(x) (results$AVE[x]-results$AVE[x-1])))
sc=results$sparsity[14:2][sc]
st=which.max(sapply(28:16,function(x) (results$AVE[x]-results$AVE[x-1])))
st=results$sparsity[28:16][st]
sm=which.max(sapply(42:30,function(x) (results$AVE[x]-results$AVE[x-1])))
sm=results$sparsity[42:30][sm]

##############################GET SELECTED FEATURES
temp=wrapper.sgcca(X=data,penalty=c(sc,st,sm,1),scale=T,
	scheme="centroid",ncomp=comps))#ncomp to explain 50% of transcripts matrix according to mfa.R
output=rbind(rowSums(do.call(rbind,temp$AVE$AVE_X)),temp$penalty)
rownames(output)=c("sum_AVE","penalty")
write.table(output,paste(subtype,BP,"descriptors",sep='.'),sep='\t',
	quote=F)
#selectVar() is the same than loadings != 0, NOT loadings.star != 0
output=lapply(temp$loadings,function(x) x[rowSums(x!=0)>0,])
lapply(1:4,function(x) 
	write.table(output[[x]],
		paste(subtype,BP,names(output)[x],sep='.'),
		sep='\t',quote=F))

#####ISSUE
#como tienes un número diferente de features para cada BP, no sabes
#si un BP tiene mejor AVE que otros, por los valores de penalty
#o porque de hecho el BP no está tan asociado con las omicas
#si fijas los penalties, igualmente las diferencias del AVE pueden
#deberse a que los penalities no son apropiados para el BP

#puedes ¿debes? repetir con los penalties ajustados y 
#matrices aleatorizadas 