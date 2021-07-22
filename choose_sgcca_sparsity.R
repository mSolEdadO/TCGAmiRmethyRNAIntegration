#aim: find a sparsity value per omic that maximizes evar &
# minimizes nfeatures
######################PIMP DATA AS NEEDED
library(data.table)
#data=fread(paste(subtype,"normalized",sep='.'))
data=fread("Her2.normalized")
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
data=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,
	function(x) t(data[x[1]:x[2],]))
names(data)=c("CpGs","transcripts","miRNAs")
tfs=read.table("Downloads/TFtargets.tsv",header=T,sep='\t')#https://github.com/slowkow/tftargets
tfs=unique(unlist(strsplit(as.character(tfs$TF),',')))
#length(tfs)
#[1] 991 lot less than https://doi.org/10.1016/j.cell.2018.01.029
#BUT tftargets does have an evidence based list of interactions
data$TFs=data$transcripts[,colnames(data$transcripts)%in%tfs]
data$transcripts=data$transcripts[,!colnames(data$transcripts)%in%tfs]
#sapply(data,dim)
#       CpGs transcripts miRNAs TFs
#[1,]     46          46     46  46
#[2,] 393132       16203    604 874

######################NARROW DOWN SPARSITY VALUES
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
#will transcripts sparsity at 1 account for corr with other blocks???
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(x,1,1,1),scale=T,scheme="centroid"))
#centroid scheme enable two components to be negatively correlated
results=list()
results$CpGs=describe(temp,grid,"CpGs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,1,x,1),scale=T,scheme="centroid"))
results$miRNAs=describe(temp,grid,"miRNAs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,1,1,x),scale=T,scheme="centroid"))
results$TFs=describe(temp,grid,"TFs")
results=do.call(rbind,results)

plots=lapply(levels(results$omic),function(i) 
	ggplot(results[results$omic==i,],aes(x=nfeatures,
		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
		scale_x_continuous(trans="log10")+
		theme(text=element_text(size=18)))
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
 #no facet_wrap coz u want indy axes
dev.off()

#turn it analyticalish by choosing the largest fall in AVE
sc=which.max(sapply(14:2,function(x) (results$AVE[x]-results$AVE[x-1])))
sc=results$sparsity[14:2][sc]
st=which.max(sapply(28:16,function(x) (results$AVE[x]-results$AVE[x-1])))
st=results$sparsity[28:16][sm]
sm=which.max(sapply(42:30,function(x) (results$AVE[x]-results$AVE[x-1])))
sm=results$sparsity[42:30][sm]
#choosen sparsity
#c(sc,st,sm)
#[1] 0.06 0.30 0.06

#######separate transcripts per enriched biological process
library(HTSanalyzeR)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
enriched=read.table("fgseaRes.tsv",sep='\t',header=T)
#enrich_GSEA.R output without leading edges
enriched=enriched[enriched$contrast=="Her2_Normal",]
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    
GS_GO_BP=GS_GO_BP[names(GS_GO_BP)%in%enriched$pathway]
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
	mart=mart)
myannot=myannot[myannot$entrezgene_id%in%unlist(GS_GO_BP),]
perBP=lapply(GS_GO_BP,function(x) 
	data$transcripts[,colnames(data$transcripts)%in%
	myannot$ensembl_gene_id[myannot$entrezgene_id%in%x]])

data$transcripts=NULL
#sapply(data,ncol)
#  CpGs miRNAs    TFs 
#393132    604   1500

##############SGCCA PER BP
test=pblapply(perBP[1:2],function(x) wrapper.sgcca(X=list(CpGs=data$CpGs,
	TFs=data$TFs,BP=x,miRNAs=data$miRNAs),penalty=c(sc,st,1,sm),
	scale=T,scheme="centroid",ncomp=15))
#ncomp to explain 50% of transcripts matrix, found with mfa.R

sapply(perBP[1:2],ncol)
#GO:0000070 GO:0000278 
#        20         90 
sapply(test,function(x) sum(x$AVE$AVE_X$BP))
#GO:0000070 GO:0000278 
# 0.9264005  0.6351158 
#does that mean GO:0000278 is not so correlated with CpG/miR/TF????
sapply(test,function(x) sapply(x$loadings,function(y) sum(rowSums(y!=0)>0)))
#       GO:0000070 GO:0000278
#CpGs        16891      16850
#TFs          1178       1160
#BP             20         90
#miRNAs        469        460
#using the same sparsity across BPs will select the same number of 
#features no mather how big is the BP

#####################CHECK LOADINGS
init=wrapper.sgcca(data,penalty=c(1,1,1),scale=T)

loadings=as.data.frame(do.call(rbind,init$loadings))
loadings$omic=substr(rownames(loadings),1,1)
loadings$omic=gsub("E","transcripts",gsub("h","miRNAs",
	gsub("c","CpGs",loadings$omic)))
loadings$omic=factor(loadings$omic,
	levels=c("CpGs","transcripts","miRNAs"))
png("loadings_l1.png")
  ggplot(loadings,aes(y=comp1,x=omic))+geom_boxplot()+
  ylab("loadings")+theme(text=element_text(size=18))
dev.off()
#transcript & miRNA loadings are normally distributed
shapiro.test(sample(loadings$comp1[loadings$omic=="transcripts"],5000))
#W = 0.99949, p-value = 0.1937
shapiro.test(loadings$comp1[loadings$omic=="miRNAs"])
#W = 0.99848, p-value = 0.8877
shapiro.test(sample(loadings$comp1[loadings$omic=="CpGs"],5000))#size must be between 3 and 5000
#W = 0.98626, p-value < 2.2e-16
###final loadgings all have p-value < 2.2e-16

################################MI NET FOR BP1
library(data.table)
library(igraph)
edges=fread("Her2.sort")
g=graph.data.frame(edges[,1:2],directed=F)
E(g)$mi=edges$V3
#drop duplicate edges
s=simplify(g,edge.attr.comb="first")
edges=as.data.frame(get.edgelist(s))
edges$MI=E(s)$mi
#drop components with no target gene
g=subgraph(s,names(compos$membership)[
	compos$membership%in%unique(compos$membership[
		names(compos$membership)%in%targes])])

i=paste(substr(edges$V1,1,1),substr(edges$V2,1,1))
methy=fread("/home/msoledad/jhu-usc.edu_BRCA.HumanMethylation450.6.lvl-3.TCGA-AO-A03N-01B-11D-A10N-05.gdc_hg38.txt")

methy_regions=apply(edges[which(i=="c c"),],1,function(x) methy$CGI_Coordinate[methy$element==x[1]]==methy$CGI_Coordinate[methy$element==x[2]])
