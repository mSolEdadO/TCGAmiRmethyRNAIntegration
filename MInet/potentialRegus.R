library(igraph)
library(data.table)
library(tidyverse)
library(ggplot2)
#sifs as graphs so nodes can easily be extracted
g=lapply(BPenriched,function(x) graph.data.frame(x[,c(1,3)],directed=F))
#########################POTENTIAL REGULATORS PER SUBTYPE#########################
nodes=lapply(g,function(x) V(x)$name)
temp=unique(unlist(nodes))
myannot=myannot[myannot$ensembl_gene_id%in%temp|myannot$mirbase_id%in%temp,]
#from https://github.com/slowkow/tftargets
tfs=fread("TFtargets.tsv")#source, tissue, TF & target hgnc symbols
#keep only those involving nodes
tfs=tfs[tfs$V3%in%myannot$hgnc_symbol&tfs$V4%in%myannot$hgnc_symbol,]
regus=sapply(nodes,function(x) table(substr(x,1,1)))
#not all transcripts are tfs,so ↓
regus[2,]=sapply(nodes,function(x) 
	sum(x%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]))
regus=data.frame(cbind(rownames(regus),regus))%>%
 pivot_longer(-V1,names_to="subtype",values_to="sum")
regus$V1=gsub("c","CpG",regus$V1)
regus$V1=gsub("h","miRNA",regus$V1)
regus$V1=gsub("E","TF",regus$V1)
regus$sum=as.numeric(as.character(regus$sum))

png("regulatorSum.png")
ggplot(regus,aes(y=sum,x=V1,fill=subtype))+
 geom_bar(position="dodge",stat="identity")+
 scale_y_continuous(trans="log10")+
 ylab("potential regulators")+
 theme(text=element_text(size=18))+xlab("")
dev.off()
#########################POTENTIAL REGULATORS PER BP#########################
#first neighbors per BP
regus=lapply(1:5,function(x) 
	lapply(GS_GO_BP[[x]],function(y) 
	 neighbors(g[[x]],which(V(g[[x]])$name%in%
	 	myannot$ensembl_gene_id[myannot$entrezgene%in%y]))$name))
#focus on potential regulators
regus=lapply(regus,function(x) 
	lapply(x,function(y) 
		y[c(which(substr(y,1,1)%in%c('c','h')),
			which(y%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]))]))
sapply(regus,function(x) sum(sapply(x,length)==0))
# Basal   Her2   LumA   LumB normal 
#    16     18      4     11     31 
#not all BPs have potential regulators

##################how many regus per bp??????????????

#data frame for barplot
names(regus)=names(top)
temp=lapply(regus,function(x) lapply(x,function(y) table(substr(y,1,1))))
temp=lapply(temp,function(x) lapply(x,function(y) cbind(names(y),y)))
temp=lapply(temp,function(x) do.call(rbind,lapply(1:length(x),function(y) 
	cbind(names(x)[y],x[[y]]))))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(temp)[x],temp[[x]]))
colnames(temp)=c("subtype","bp","omic","sum")
temp$sum=as.numeric(as.character(temp$sum))

ggplot(temp,aes(x=subtype,y=x,fill=regulator))+
geom_bar(position="stack", stat="identity")+
ylab("biological processes")+xlab("")+
theme(text=element_text(size=18))+
scale_fill_manual(values=c("skyblue4","skyblue3","skyblue2","skyblue1","plum3","plum4"))
#colors have to be adjusted
#for heatmap
temp=lapply(regus,function(x) lapply(x,function(y) table(substr(y,1,1))))
temp=lapply(temp,function(x) lapply(x,function(y) cbind(names(y),y)))
temp=lapply(temp,function(x) do.call(rbind,lapply(1:length(x),function(y) 
	cbind(names(x)[y],x[[y]]))))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(temp)[x],temp[[x]]))
colnames(temp)=c("subtype","bp","omic","sum")
temp$sum=as.numeric(as.character(temp$sum))
M=graph.data.frame(temp[,1:2],directed=F)
E(M)$weight=temp$sum
M=as.matrix(M[unique(temp[,2]),unique(temp[,1])])
M[M==0]=NA
M=lapply(names(top),function(x) M[,grep(x,colnames(M))])
total=sapply(M,rowSums)
M=do.call(cbind,lapply(1:5,function(x) 100*M[[x]]/total[,x]))
#u still have get % of to plot
mir=fread("miR.ids.map.tsv",skip=1)