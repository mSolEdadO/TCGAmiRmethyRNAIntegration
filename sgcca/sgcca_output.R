#######################net at different cor.value
library(data.table)
library(igraph)
library(biomaRt)
library(tidyverse)

files=list.files()
files=files[grep("net",files)]
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
descri=lapply(files,function(y) {
	net=fread(y);
	BP=gsub(".net","",gsub("Her2.","",y))
	geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),
		filters = 'go',values = BP, mart = mart);
	geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP]);
	size=sapply(seq(0,.9,.1),function(x) {
	g=graph.data.frame(net[net$corr>x,1:2],directed=F);
	edges=ecount(g);
	nodes=V(g)$name;
	bpNodes=sum(nodes%in%geneSet)
	descri=list(edges=edges,CpGs=sum(substr(nodes,1,1)=="c"),
		TFs=sum(substr(nodes,1,1)=="E")-bpNodes,
		miRNAs=sum(substr(nodes,1,1)=="h"),
		BP=bpNodes)
	return(descri)})
})
descri=do.call(rbind,lapply(1:length(descri),function(x) 
	cbind(gsub(".net","",gsub("Her2.","",files[x])),
		rownames(descri[[x]]),descri[[x]])))
colnames(descri)=c("BP","feature",seq(0,0.9,0.1))
rownames(descri)=NULL
descri=as.data.frame(as.matrix(descri))
descri=as.data.frame(descri%>%pivot_longer(-c(1,2),names_to="cor",
	values_to="nfeature"))

####################choose BPs to describe
library(ggplot2)

files=list.files()
files=files[grep("descri",files)]
data=lapply(files,read.table)
#verify is the same order
#sum(sapply(strsplit(files,".",fixed=T),function(x) x[2])==unique(descri$BP))
#df to plot
temp=sapply(unique(descri$BP),function(x) 
	sum(descri$nfeature[descri$BP==x&descri$feature!="edges"&descri$cor==0]))
temp=as.data.frame(cbind(sapply(data,function(x) x[1,4]),temp))
colnames(temp)=c("sum_AVE","nfeatures")
png("quadrants.png")
ggplot(temp,aes(x=nfeatures,y=sum_AVE))+geom_vline(aes(xintercept=1e4))+
geom_hline(aes(yintercept=.5))+geom_point()+
scale_x_continuous(trans="log10")+theme(text=element_text(size=18))
dev.off()
#most BPs are fairly explained by selected features
#since ncomp is the same, non-explained BPs may depend on other omics
cor.test(temp$sum_AVE,temp$nfeatures,method="spearman")#best method???
#S = 22761, p-value = 0.03051
#       rho 
#-0.3159477 
#small but significant cor between axes

################### choose edges?
lala=descri[descri$BP%in%levels(descri$BP)[1:6],]
##I want only the most important correlations BUT, throwing edges
# also throws BP transcripts
png("Cov.vs.Edges.png")
 ggplot(lala,aes(x=cor,y=nfeature,color=feature))+geom_line()+
 scale_y_continuous(trans="log10")+geom_point()+facet_wrap(~BP,ncol=3,
 nrow=2)+theme(text=element_text(size=18))
dev.off()
#so I CAN NOT choose edges based on correlation
