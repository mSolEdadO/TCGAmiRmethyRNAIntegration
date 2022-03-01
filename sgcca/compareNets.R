#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
fun=commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(igraph)
library(ape)
library(ggplot2)
#########################LOAD NETS
files=list.files()
files=files[grep(fun,files)]
files=files[grep("tsv",files)]
edges=lapply(files,read_tsv,col_names=F)
names(edges)=sapply(strsplit(files,".",fixed=T),function(x) x[1])

nets=lapply(edges,function(x) {g=graph.data.frame(x[,1:2]);#directed graph
							   E(g)$weight=x$X4;
							   return(g)})
#########################WHICH SUBTYPES ARE MORE ALIKE?
d=sapply(nets,function(x) 
	sapply(nets,function(y) 
		ecount(intersection(x,y))/ecount(graph.union(x,y))))
tree=hclust(as.dist(1-d))
#to build a consensus with all the functions
write.tree(as.phylo(tree),paste(fun,"newick",sep='.'))
#if u wanna plot
#png(paste(fun,"tile.png",sep='.'))
#cbind(rownames(d),d)%>%as.data.frame%>%
#	pivot_longer(-1,names_to="subtype",values_to="ecount")%>%
#	ggplot(aes(V1,subtype,fill=log10(as.numeric(ecount))))+
#	geom_tile()+scale_fill_viridis_c(name="log(intersection)")+
#	xlab("")+ylab("")+theme(text=element_text(size=18))
#dev.off()

#########################WHICH TYPE OF EDGES ARE MORE SHARED
edges=do.call(rbind,lapply(1:5,function(x)
	cbind(subtype=names(edges)[x],edges[[x]])))
#get edge type
temp=t(apply(edges,1,function(x)
	cbind(substr(x[2],1,1),substr(x[3],1,1))))
edges$type=apply(temp,1,function(y)
	#so E-h is the same as h-E
	paste(sort(y),collapse='-'))
edges$i=paste(edges$X1,edges$X2)
shared=edges%>%group_by(type)%>%count(i)%>%group_map(~table(.x$n>1))
names(shared)=edges%>%group_by(type)%>%group_keys()%>%unlist()
ps=p.adjust(sapply(names(shared),function(z) {#per edge type 
	mat=rbind(shared[[z]][c("FALSE","TRUE")],#bind even when T/F vals may be missing, with the
		sapply(c("FALSE","TRUE"),function(y)#sum of all exclusive/shared
		 sum(sapply(shared[names(shared)!=z],#of the other types
		 	function(x) x[y]),na.rm=T)));
	colnames(mat)=c("FALSE","TRUE")#correct NA
	mat[is.na(mat)]=0#correct NA
	fisher.test(mat,alternative='l')$p.val}))#test
print(paste(names(which(ps<0.01)),
	"edges are more shared than expected, with p.values",
	ps[ps<0.01],sep=' '))
#alternative
#temp=chisq.test(do.call(rbind,shared))
#print(temp$residuals)
#but this way u dont get a p.value per edge type

#########################INTERSECT NETS
shared=intersection(nets$Basal,
					nets$Her2,
					nets$LumA,
					nets$LumB,
					nets$Normal)
n=5
if(ecount(shared)==0){
	shared=intersection(nets$Basal,
						nets$Her2,
						nets$LumA,
						nets$LumB)
	n=4
}
shared=induced_subgraph(shared,which(degree(shared)>0))#drop loose nodes

#DO WEIGHTS CHANGE ACROSS SUBTYPES?
ws=edge_attr(shared)
ws=as.data.frame(do.call(cbind,ws[1:n]))
colnames(ws)=names(nets)[1:n]
ws$edge=1:nrow(ws)
#png(paste(fun,"weights.png",sep='.'))
ws%>%pivot_longer(-edge,names_to="subtype",values_to="weight")%>%
	ggplot(aes(subtype,edge,fill=weight))+geom_tile()+
	scale_fill_viridis_c()+ylab("")+xlab("")+
	theme(text=element_text(size=18),axis.ticks.y=element_blank(),
		axis.text.y=element_blank())
#dev.off()

#########################ADD KNOWN REGULATORY EDGES
files=list.files()
files=files[grep(fun,files)]
files=files[grep("moti",files)]
known=lapply(files,read_tsv,col_names=F)
known=do.call(rbind,known)
known=graph.data.frame(known[,1:2])
#focus only on the nodes with shared links 
known=make_ego_graph(known,
	nodes=which(V(known)$name%in%V(shared)$name))
final=rbind(get.edgelist(shared),
			unique(do.call(rbind,lapply(known,get.edgelist))))
nrow(final)
[1] 206
#por gene names y agrupa targets por función?????????


temp=graph.edgelist(final)
plot(temp,vertex.color="gray",vertex.label.cex=1.5,
	vertex.frame.color=NA,vertex.label.family="sans",
	vertex.label.color="black")
#WHAT MAKES SPECIAL THIS EDGES?????????

