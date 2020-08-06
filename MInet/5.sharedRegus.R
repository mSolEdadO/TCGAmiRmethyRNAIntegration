library(ggplot2)
#######################REGULATORS SHARED WITHIN SUTYPES#########################
#regulators shared by BP
temp=lapply(regus,function(x) table(unlist(x)))
#data frame to plot
temp=lapply(temp,function(x) cbind(substr(names(x),1,1),x))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
 cbind(names(temp)[x],temp[[x]]))))
temp$V2=gsub("c","CpG",temp$V2)
temp$V2=gsub("h","miRNA",temp$V2)
temp$V2=gsub("E","TF",temp$V2)
colnames(temp)=c("subtype","regulator","bp")
temp$bp=as.numeric(as.character(temp$bp))
temp$regulator=factor(temp$regulator,levels=c("CpG","TF","miRNA"))

png("BPpeRegu.png")
 ggplot(temp,aes(y=bp,x=subtype,color=subtype))+geom_boxplot()+
 facet_wrap(~regulator)+ylim(0,20)+xlab("")+ylab("processes by regulator")+
 theme(text=element_text(size=18),axis.text.x=element_blank())
dev.off()

#exclusivity of regulators
#get categories needed
temp$exclu=temp$bp==1
temp$exclu=gsub(TRUE,"exclusive",temp$exclu)
temp$exclu=gsub(FALSE,"shared",temp$exclu)

png("regulatorShared.png")
 ggplot(temp,aes(y=bp,fill=regulator,x=subtype))+
 geom_bar(position="fill", stat="identity")+
 facet_wrap(~exclu)+xlab("")+ylab("proportion")+
 theme(text=element_text(size=18),axis.text.x=element_text(angle=45))+
 scale_fill_manual(values=gray.colors(5))
dev.off()

#fisher test per omic
temp=lapply(names(regus),function(x) table(temp[temp$subtype==x,c(2,4)]))
names(temp)=names(regus)
temp=lapply(c("CpG","TF","miRNA"),function(x) 
	lapply(temp,function(y) rbind(y[rownames(y)==x,],colSums(y[rownames(y)!=x,]))))
names(temp)=c("CpG","TF","miRNA")
sapply(temp,function(x) sapply(x,function(y) p.adjust(fisher.test(y,alternative="g")$p.val,"fdr")))
#                          CpG          TF      miRNA
#Basal.exclusive  9.813473e-01 0.026251494 0.58243408
#Her2.exclusive   9.997335e-01 0.001561389 0.47936372
#LumA.exclusive   9.894221e-01 0.746753451 0.09308053
#LumB.exclusive   8.887973e-15 1.000000000 0.99279081
#normal.exclusive 5.617978e-01 0.989768640 0.03146967
sapply(temp,function(x) sapply(x,function(y) p.adjust(fisher.test(y,alternative="l")$p.val,"fdr")))
#                          CpG           TF      miRNA
#Basal.exclusive  0.0372178873 9.882466e-01 0.56898844
#Her2.exclusive   0.0008306361 9.994873e-01 0.66241785
#LumA.exclusive   0.0631914348 4.093731e-01 0.96041890
#LumB.exclusive   1.0000000000 1.490616e-11 0.01802077
#normal.exclusive 1.0000000000 2.310520e-02 0.98554771


#% of regulators shared by BP
temp=lapply(regus,function(x) x[sapply(x,length)>0])
temp=lapply(temp,function(x) do.call(rbind,lapply(1:length(x),function(y) 
 cbind(names(x)[y],x[[y]]))))
#matrices per regulator type
temp=lapply(temp,function(x) 
 as.matrix(graph.edgelist(x,directed=F)[unique(x[,1]),unique(x[,2])])) 
temp=lapply(temp,function(x) lapply(c("c","E","h"),function(y)
 x[,substr(colnames(x),1,1)%in%y]))
#matrices of totals shared
temp=lapply(temp,function(x) lapply(x,function(y) y%*%t(y)))
#percentage
temp=lapply(temp,function(x) lapply(x,function(y) y/diag(y)))
#data frame to plot
temp=lapply(temp,function(x) lapply(x,function(y) y[upper.tri(y)]))
temp=lapply(temp,function(x) do.call(rbind,lapply(1:3,function(y) 
	cbind(c("CpG","TF","miRNA")[y],x[[y]]))))
withinSubty=data.frame(do.call(rbind,lapply(1:5,function(x)
 cbind(names(temp)[x],temp[[x]]))))
#include no regulators = genes annotated in the functions
temp1=lapply(GS_GO_BP,function(x) lapply(x,function(y) 
	myannot$ensembl_gene_id[myannot$entrezgene%in%y]))
temp1=lapply(1:5,function(x) lapply(temp1[[x]],function(y) y[y%in%V(g[[x]])$name]))
temp1=lapply(1:5,function(x) lapply(temp1[[x]],function(y)
	y[!y%in%unlist(regus)]))
#repeat % estimation
temp1=lapply(temp1,function(x) do.call(rbind,lapply(1:length(x),function(y) 
 cbind(names(x)[y],x[[y]]))))
temp1=lapply(temp1,function(x) 
 as.matrix(graph.edgelist(x,directed=F)[unique(x[,1]),unique(x[,2])])) 
temp1=lapply(temp1,function(x) x%*%t(x))
temp1=lapply(temp1,function(x) x/diag(x))
temp1=lapply(temp1,function(x) x[upper.tri(x)])
names(temp1)=names(regus)
temp1=data.frame(do.call(rbind,lapply(1:5,function(x)
 cbind(names(temp1)[x],"transcript",temp1[[x]]))))
#paste all together
withinSubty=data.frame(rbind(withinSubty,temp1))
colnames(withinSubty)=c("subtype","omic","shared")
withinSubty$shared=100*as.numeric(as.character(withinSubty$shared))

png("perSubtypeSharing.png")
 ggplot(temp,aes(x=shared))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~omic,ncol=1)+xlab("% shared")+theme(text=element_text(size=18))
dev.off()

#######################REGULATORS SHARED BETWEEN SUTYPES#########################
#keep only BP with associated regulators
temp=lapply(regus,function(x) x[sapply(x,length)>0])
#get BPs enriched in several nets
i=names(which(table(unlist(lapply(temp,names)))>1))
temp=lapply(i,function(x) sapply(regus,function(y) y[names(y)==x]))
temp=lapply(temp,function(y) y[sapply(y,function(x) length(unlist(x)))>0])
#make an edgelist of subtype vs regulator per BP
temp=lapply(1:length(i),function(y) 
	do.call(rbind,lapply(1:length(temp[[y]]),function(x) 
		cbind(names(temp[[y]])[x],temp[[y]][[x]][[1]]))))
names(temp)=i
#get matrices of subtype vs regulator per BP
temp=lapply(temp,function(x) 
	as.matrix(graph.edgelist(x,directed=F)[unique(x[,1]),unique(x[,2])]))

#how many BPs are shared between subtypes
temp1=table(sapply(temp,function(x) 
	paste(sapply(strsplit(rownames(x),".",fixed=T),function(y) y[1]),collapse='+')))
#data frame to plot
temp1=data.frame(cbind(names(temp1),temp1))
temp1$temp1=as.numeric(as.character(temp1$temp1))
temp1=temp1[order(temp1$temp1,decreasing=T),]
temp1$V1=factor(temp1$V1,levels=as.character(temp1$V1))

png("BPonSubt.png")
 ggplot(temp1,aes(x=temp1,y=V1))+geom_bar(position="dodge", stat="identity")+
 ylab("")+xlab("count")+theme(text=element_text(size=18))
dev.off()

#matrix per omic for each BP
temp=lapply(temp,function(x) lapply(c("c","E","h"),function(y) x[,substr(colnames(x),1,1)==y]))
temp=lapply(temp,function(x) lapply(x,function(y) y%*%t(y)))
temp=lapply(temp,function(x) lapply(x,function(y) y/diag(y)))
#data frame the matrices
temp=lapply(1:length(temp),function(z) lapply(1:3,function(y) 
	do.call(rbind,lapply(1:ncol(temp[[z]][[1]]),function(x) 
		cbind(colnames(temp[[z]][[y]])[x],temp[[z]][[y]][x,-x])))))
temp=lapply(1:3,function(x) do.call(rbind,lapply(temp1,function(y) y[[x]])))
#data frame to plot
acrosSubty=data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(c("CpG","TF","miRNA")[x],temp[[x]]))))
colnames(acrosSubty)=c("omic","subtype","shared")
acrosSubty$subtype=sapply(strsplit(as.character(acrosSubty$subtype),".",fixed=T),
	function(x) x[1])
acrosSubty$shared=100*as.numeric(as.character(acrosSubty$shared))
acrosSubty=acrosSubty[!is.na(acrosSubty$shared),]
acrosSubty$omic=factor(acrosSubty$omic,levels=c("CpG","TF","miRNA"))

png("accrosSubtypeSharing.png")
 ggplot(acrosSubty,aes(x=shared))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~omic,ncol=1)+theme(text=element_text(size=18))
dev.off()