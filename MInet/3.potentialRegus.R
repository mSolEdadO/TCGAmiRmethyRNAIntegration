library(igraph)
library(data.table)
library(tidyverse)
library(ggplot2)
#sifs as graphs so nodes can easily be extracted
g=lapply(BPenriched,function(x) graph.data.frame(x[,c(1,3)],directed=F))
#########################POTENTIAL REGULATORS PER SUBTYPE#########################
nodes=lapply(g,function(x) V(x)$name)
#from https://github.com/slowkow/tftargets
tfs=fread("TFtargets.tsv")#source, tissue, TF & target hgnc symbols
regus=sapply(nodes,function(x) table(substr(x,1,1)))
#not all transcripts are tfs,so ↓
regus[2,]=sapply(nodes,function(x) 
	sum(x%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]))
regus=data.frame(cbind(rownames(regus),regus))
regus=regus%>%pivot_longer(-V1,names_to="subtype",values_to="sum")
regus$V1=gsub("c","CpG",regus$V1)
regus$V1=gsub("h","miRNA",regus$V1)
regus$V1=gsub("E","TF",regus$V1)
regus$sum=as.numeric(as.character(regus$sum))
regus$V1=factor(regus$V1,levels=c("CpG","TF","miRNA"))

png("regulatorSum.png")
ggplot(regus,aes(y=sum,x=V1,fill=subtype))+
 geom_bar(position="dodge",stat="identity")+
 scale_y_continuous(trans="log10")+
 ylab("potential regulators")+
 theme(text=element_text(size=18))+xlab("")
dev.off()
#########################POTENTIAL REGULATORS PER BP#########################
#first neighbors per BP
regus=lapply(1:5,function(z) 
	lapply(GS_GO_BP[[z]],function(y) 
		unique(unlist(sapply(myannot$ensembl_gene_id[myannot$entrezgene%in%y],
			function(x) neighbors(g[[z]],which(V(g[[z]])$name==x))$name)))))
#focus on potential regulators
regus=lapply(regus,function(x) 
	lapply(x,function(y) 
		y[c(which(substr(y,1,1)%in%c('c','h')),
			which(y%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]))]))
names(regus)=names(top)
#not all BPs have potential regulators
sapply(regus,function(x) sum(sapply(x,length)==0))
# Basal   Her2   LumA   LumB normal 
#     0      0      0      0      0 
#not interested in BPs enriched only in the net from normal tissue
regus$normal=regus$normal[names(regus$normal)%in%unlist(sapply(regus[1:4],names))]

#regulator's combination per subtype
temp=lapply(regus,function(x) table(unlist(lapply(x,function(y) 
	paste(unique(substr(y,1,1)),collapse="")))))
temp=lapply(temp,function(x) cbind(names(x),x))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(temp)[x],temp[[x]]))))
colnames(temp)=c("subtype","combo","sum")
temp$sum=as.numeric(as.character(temp$sum))
#make omic ids explicit
temp$combo=gsub("c","CpG+",temp$combo)
temp$combo=gsub("h","miRNA+",temp$combo)
temp$combo=gsub("E","TF",temp$combo)
temp$combo[temp$combo=="miRNA+CpG+"]="CpG+miRNA"
temp$combo[temp$combo=="miRNA+CpG+TF"]="CpG+TF+miRNA"
temp$combo[temp$combo=="CpG+miRNA+TF"]="CpG+TF+miRNA"

png("regulatorsCombo.png")
 ggplot(temp,aes(x=subtype,y=sum,fill=combo))+
 geom_bar(position="stack", stat="identity")+
 ylab("biological processes")+xlab("")+
 theme(text=element_text(size=18))+
 scale_fill_manual("regulators",values=gray.colors(5))
dev.off()
#colors have to be adjusted

#% of each potential regulator per BP
temp=lapply(regus,function(x) lapply(x,function(y) table(substr(y,1,1))))
temp=lapply(temp,function(x) lapply(x,function(y) cbind(names(y),y)))
temp=lapply(temp,function(x) do.call(rbind,lapply(1:length(x),function(y) 
	cbind(names(x)[y],x[[y]]))))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(temp)[x],temp[[x]]))))
colnames(temp)=c("subtype","bp","omic","sum")
temp$sum=as.numeric(as.character(temp$sum))
temp$subtype=paste(temp$subtype,temp$omic)
#matrix of subtype-regulator type vs BP
M=graph.data.frame(temp[,1:2],directed=F)
E(M)$weight=temp$sum
M=as.matrix(M[unique(temp[,2]),unique(temp[,1])])
M=lapply(names(top),function(x) M[,grep(x,colnames(M))])
total=sapply(M,rowSums,na.rm=T)
#get percentages
M=do.call(cbind,lapply(1:5,function(x) 100*M[[x]]/total[,x]))
#order columns
M=M[,c(1:12,15,13,14)]
M[is.na(M)]=0

pdf("RegusTot.pdf",height=18)
gplots::heatmap.2(M,Colv=F,dendrogram='r',scale='n',sepcolor="gray",
	col=rev(heat.colors(100)),trace='n',labCol=rep(c("CpG","TF","miRNA"),5),
	lmat = rbind(c(0,3),c(2,1),c(0,4)),lhei=c(0.2,5,2),lwid=c(1,6),
	key.title="",margins=c(2,25),density.info="n",
	add.expr=text(x=c(2,5,8,11,14),y=178,xpd=NA,labels=names(top)),key.xlab="%")
#xy positions need to be adjusted manually
dev.off()

#regulators per BP
temp$subtype=sapply(strsplit(temp$subtype," "),function(x) x[1])
temp$omic=gsub("c","CpG",temp$omic)
temp$omic=gsub("h","miRNA",temp$omic)
temp$omic=gsub("E","TF",temp$omic)
temp$omic=factor(temp$omic,levels=c("CpG","TF","miRNA"))

png("perBPdistris.png")
 ggplot(temp,aes(y=sum,x=subtype,color=subtype))+ylim(0,900)+geom_boxplot()+
 facet_wrap(~omic)+ylab("regulators per process")+
 theme(text=element_text(size=18),axis.text.x=element_blank())+xlab("")
dev.off()
#Removed 1 rows containing non-finite values (stat_boxplot). 
#test distributions
temp=lapply(unique(temp$omic),function(x) temp[temp$omic==x,])
names(temp)=sapply(temp,function(x) x$omic[1])
lapply(temp,function(z) 
	matrix(round(p.adjust(sapply(names(top),function(x) sapply(names(top),
		function(y) ks.test(z$sum[z$subtype==x],z$sum[z$subtype==y])$p.val)),
					"fdr"),4),ncol=5))
#mir=fread("miR.ids.map.tsv",skip=1)