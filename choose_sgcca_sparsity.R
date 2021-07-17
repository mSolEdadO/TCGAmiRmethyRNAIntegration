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
tfs=readLines("TFs_Ensembl_v_1.01.txt")
#https://www.cell.com/cell/fulltext/S0092-8674(18)30106-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418301065%3Fshowall%3Dtrue
data$TFs=data$transcripts[,colnames(data$transcripts)%in%tfs]

#separate transcripts per enriched biological process
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
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(x,1,1),scale=T,scheme="centroid"))
results=list()
results$CpGs=describe(temp,grid,"CpGs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,x,1),scale=T,scheme="centroid"))
results$miRNAs=describe(temp,grid,"miRNAs")
temp=pblapply(grid,function(x) 
	wrapper.sgcca(data,penalty=c(1,1,x),scale=T,scheme="centroid"))
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
results$sparsity[14:2][s]
#[1] 0.04"
#choosen sparsity
sm=which.max(sapply(28:16,function(x) (results$AVE[x]-results$AVE[x-1])))
results$sparsity[28:16][s]
#[1] 0.2
st=which.max(sapply(42:30,function(x) (results$AVE[x]-results$AVE[x-1])))
results$sparsity[42:30][s]
#[1] 0.2

##############SGCCA PER BP
test=pblapply(perBP[1:2],function(x) wrapper.sgcca(X=list(CpGs=data$CpGs,
	TFs=data$TFs,BP=x,miRNAs=data$miRNAs),penalty=c(sc,st,1,sm),
	scale=T,scheme="centroid",ncomp=15))
#ncomp to explain 50% of transcripts matrix, found with mfa.R

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

################################COMPARE INIT & FINAL
#%of selected features
100*sapply(final$loadings,function(x) 
	sum(x!=0))/sapply(init$loadings,function(x) sum(x!=0))
#       CpGs transcripts      miRNAs 
#   6.187235   63.565029   46.026490 
#% of AVE from initial
100*do.call(rbind,final$AVE[[1]])/do.call(rbind,init$AVE[[1]])
#                comp1
#CpGs         80.10103
#transcripts  96.51604
#miRNAs      103.08098
#####whats with miRNAs?????
temp=as.data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(init$loadings[[x]],final$loadings[[x]]))))
colnames(temp)=c("initial","final")
temp$omic=substr(rownames(temp),1,1)
temp$omic=gsub("E","transcripts",gsub("h","miRNAs",
 gsub("c","CpGs",loadings$omic)))
plots=lapply(unique(temp$omic),function(x) 
	ggplot(temp[temp$omic==x,],aes(x=initial,y=final))+geom_point()+
	theme(text=element_text(size=18)))
png("loadings_change.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

