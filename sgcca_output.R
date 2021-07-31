library(ggplot2)

####################choose BPs to describe
files=list.files()
files=files[grep("descri",files)]
data=lapply(files,read.table)
#for x in $(ls|grep 'net');do echo $x; cut -f1,2 $x|perl -pe 's/\t/\n/'|sort|uniq|wc -l;done
size=read.table("sizes")
#verify is the same order
#sum(sapply(strsplit(files,".",fixed=T),function(x) x[2])==size$V1)
#df to plot
temp=as.data.frame(cbind(sapply(data,function(x) x[1,4]),size$V2))
colnames(temp)=c("sum_AVE","nfeatures")
png("quadrants.png")
ggplot(temp,aes(x=nfeatures,y=sum_AVE))+geom_vline(aes(xintercept=1e4))+
geom_hline(aes(yintercept=.5))+geom_point()+
scale_x_continuous(trans="log10")+theme(text=element_text(size=18))
dev.off()
#most BPs are fairly explained by selected features
#since ncomp is the same, non-explained BPs may depend on other omics
cor.test(temp$sum_AVE,temp$nfeatures,method="spearman")#best method???
#S = 22636, p-value = 0.03471
#       rho 
#-0.3087687 
#small but significant cor between axes

#######################net at different cor.value
library(data.table)
library(igraph)
library(biomaRt)

files=list.files()
files=files[grep("net",files)]
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
lapply(files,function(y) {
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