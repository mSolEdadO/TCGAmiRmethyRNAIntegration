#if input is from cluster
#files=list.files()
#files=files[grep("params",files)]
#temp=lapply(files,read.table,head=T)
#temp=do.call(rbind,temp)

#if input is from your lap
temp=read.table("penalty_search.tsv",sep='\t',header=T)

#######################################DIAGNOSTIC PLOTS
library(ggplot2)

temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
png("AVE.png",width=800)
 ggplot(temp,aes(y=AVE,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))
dev.off()
png("nfeatures.png",width=800)
 ggplot(temp,aes(y=nfeatures,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))+
 scale_y_continuous(trans="log10")
dev.off()
#plot median AVE vs meadian nfeatures
omics=levels(temp$omic)
omics=lapply(omics,function(x) temp[temp$omic==x,])
omics=lapply(omics,function(x) as.data.frame(apply(x,2,as.numeric)))
names(omics)=levels(temp$omic)
omics=lapply(omics,function(y) sapply(unique(y$penalty),function(x) 
apply(y[y$penalty==x,],2,median)))
omics=lapply(omics,function(x) as.data.frame(t(x)))
#indi plots or CpGs will determine axis
plots=lapply(1:3,function(x) ggplot(omics[[x]],
	aes(x=nfeatures,y=AVE,col=penalty))+geom_point()+
	ggtitle(names(omics)[x])+theme(text=element_text(size=18))+
	scale_x_continuous(trans="log10")+geom_line())
png("sparsity_search_0.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

#######################################PENALTIES SUGGESTED BY PLOTS
grid=seq(0.01,0.1,0.01)
#grid=seq(0.01,0.9,length=10)# first grid, tested in the cluster
slopes=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
	(y$AVE[x+1]-y$AVE[x])/(y$nfeatures[x+1]-y$nfeatures[x])))
slopes$miRNAs[abs(slopes$miRNAs)=="Inf"]=NA
penalty=sapply(slopes,function(x) grid[which.max(x)+1])
#the first grid selects always 0.1088889 as penalty
#so it was necesary to check for a grid with smaller values
#2nd grid selects
#       CpGs transcripts      miRNAs 
#       0.05        0.03        0.08 

#######################################FINAL SGCCA
library(igraph)
library(mixOmics)
library(data.table)

data=fread(paste(subtype,"normalized",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
n=ncol(data)
subn=round(n*.5)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")

final=wrapper.sgcca(X=data,penalty=penalty,scale=T,
	scheme="centroid",ncomp=15)#ncomp to explain 50% of transcripts matrix according to mfa.R
output=rbind(rowSums(do.call(rbind,temp$AVE$AVE_X)),temp$penalty)
rownames(output)=c("sum_AVE","penalty")
#             CpGs transcripts    miRNAs
#sum_AVE 0.4714325   0.4188032 0.4021734
#penalty 0.0500000   0.0300000 0.0800000

#######################################CHECK ENRICHMENTS
library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(HTSanalyzeR)
library("stringr")
library(GO.db)

#get transcripts per component
selected=apply(final$loadings$transcripts,2,function(x) 
	rownames(final$loadings$transcripts)[x!=0])
#get dbs 
GS_KEGG<-KeggGeneSets(species = "Hs")
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    
diccionario_kegg<-as.list(KEGGPATHID2NAME)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol",
	"entrezgene_id"), mart=mart)
universo=as.character(myannot$entrezgene_id[!is.na(myannot$entrezgene_id)])
#transform to entrez id
sets=lapply(selected,function(y) myannot$entrezgene_id[myannot$ensembl_gene_id%in%y])
sets=lapply(sets,function(x) x[!is.na(x)])
enrichBP=lapply(sets,function(y)
	as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_BP,
 	universe= universo, 
 	hits=as.character(y),
 	minGeneSetSize = 15, 
 	pAdjustMethod = "fdr"))
sapply(enrichBP,function(x) sum(x$Adjusted.Pvalue<0.05))
#no positive results
#observed hits per component stay<6
sapply(enrichKEGG,function(x) max(x[,5]))
 comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 comp11 
     3      2      2      4      3      2      3      2      4      2      1 
comp12 comp13 comp14 comp15 
     2      3      2      3 
> sapply(enrichBP,function(x) max(x[,5]))
 comp1  comp2  comp3  comp4  comp5  comp6  comp7  comp8  comp9 comp10 comp11 
     3      2      2      2      4      2      2      3      5      2      2 
comp12 comp13 comp14 comp15 
     3      2      2      3 