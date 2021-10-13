#for the cluster
#calculates the correlation between the ensembl id in args[2] 
#and all the rows after it in the matrix

#!/usr/bin/env Rscript

########BUILDED A MATRIX PER CLUSTER WITH THIS
#gps=read.csv("4clusters.csv")
#data=fread("../Basal.mtrx")
#gps=read.csv("4clusters.csv")
#data=data[grep("ENSG",data$V1),]
#data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#data=lapply(unique(gps$group),function(x) 
#	data[,colnames(data)%in%gps$sample_id[gps$group==x]])
#names(data)=c("clust2","clust3","clust1","clust4")
#lapply(1:4,function(x) write.table(data[[x]],
#	paste(names(data)[x],"mtrx",sep='.'),sep='\t',quote=F))

library(data.table)

args=commandArgs(trailingOnly=TRUE)
gpo=args[1]#name of a cluster matriX
gene=args[2]#ENSG id to test
data=fread(paste(gpo,"mtrx",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
i=which(rownames(data)==gene)

#the last gene doesn't need any calculation, it already has been paired
edges=as.data.frame(t(sapply((i+1):nrow(data),function(x) 
 cbind(rownames(data)[c(i,x)],cor.test(data[i,],data[x,])[3:4]))))
colnames(edges)=c("node1","node2","pval","rho")
edges$adjustedPval=p.adjust(edges$pval,"fdr")
#wont p.adjust fall short when joining all files?
edges=apply(edges,2,unlist)
write.table(edges[edges$adjustedPval<0.01,],paste(gpo,gene,
	"edges",sep='.'),sep='\t',quote=F,row.names=F)
