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
if(sum(ps<0.01)>0){
print(paste(names(which(ps<0.01)),
	"edges are more shared than expected, with p.values",
	ps[ps<0.01],sep=' '))
}
#alternative
#temp=chisq.test(do.call(rbind,shared))
#print(temp$residuals)
#but this way u dont get a p.value per edge type

#########################EXCLUSIVE EDGES
exclusive=edges[!edges$i%in%edges$i[duplicated(edges$i)],]
#net with exclusive edges
g_exclusive=exclusive%>%group_by(subtype)%>%
			group_map(~graph.data.frame(.x[,c("X1","X2")]))
names(g_exclusive)=exclusive%>%group_by(subtype)%>%group_keys()%>%unlist()
dg=sapply(g_exclusive,degree)
#df to plot
temp=as.data.frame(do.call(rbind,sapply(1:5,function(x) 
	cbind(subtype=names(dg)[x],g=dg[[x]]))))
png(paste(fun,"degree.png",sep='_'))
ggplot(temp,aes(x=as.numeric(g),col=subtype,fill=subtype))+
	geom_density(aes(y=..scaled..),alpha=0.3)+
	scale_x_continuous(trans="log10")+xlab("degree")+
	theme(text=element_text(size=18))
dev.off()
#check high degree nodes
temp=sapply(dg,function(x) names(which(x>mean(x))))#1:present,0:absent
temp=table(as.data.frame(do.call(rbind,sapply(1:5,function(x)
 cbind(temp[[x]],names(temp)[x])))))
temp=apply(temp,1,paste,collapse="")#0,1 chains ordered as the nets
if(sum(temp=="11110")){
	print(paste(names(which(temp=="11110")),
	"is a high degree node within exclusive edges of the 4 subtype networks"))
}
if(sum(temp=="11111")){
	print(paste(names(which(temp=="11111")),
	"is a high degree node within exclusive edges of the 5 networks"))
}


#########################INTERSECT NETS
n=5
shared=intersection(nets$Basal,
					nets$Her2,
					nets$LumA,
					nets$LumB,
					nets$Normal)
if(ecount(shared)==0){
shared=intersection(nets$Basal,
					nets$Her2,
					nets$LumA,
					nets$LumB)
	if(ecount(shared)==0){
		stop("No shared edges")
	}
n=4
print("Only the networks of the subtypes share any edges")
}
shared=induced_subgraph(shared,which(degree(shared)>0))#drop loose nodes

#DO WEIGHTS CHANGE ACROSS SUBTYPES?
ws=edge_attr(shared)
ws=as.data.frame(do.call(cbind,ws[1:n]))
colnames(ws)=names(nets)[1:n]
#ws$edge=1:nrow(ws)
#png(paste(fun,"weights.png",sep='.'))
#ws%>%pivot_longer(-edge,names_to="subtype",values_to="weight")%>%
#	ggplot(aes(subtype,edge,fill=weight))+geom_tile()+
#	scale_fill_viridis_c()+ylab("")+xlab("")+
#	theme(text=element_text(size=18),axis.ticks.y=element_blank(),
#		axis.text.y=element_blank())
#dev.off()

#plot
demir=read_tsv("DE.miR.tsv")
de=read_tsv("DE.genes.tsv")
dm=read_tsv("DMcpgs-RUV.tsv")
da=list(cpgs=dm,genes=de,mir=demir)
colnames(da$genes)[1]="id"

da=lapply(da,function(x) x[which(x$id%in%V(shared)$name),])
da=as.data.frame(unique(mapply(rbind,da$cpgs[,c(1,2,7)],
	do.call(rbind,lapply(da[2:3],function(x) x[,c(1,2,5)])))))
da=da[order(match(da$id,V(shared)$name)),]
colnames(da)[3]="fc"
lay=layout.auto(shared)
cols=rev(RColorBrewer::brewer.pal(name="RdBu",n=10))
da$cols=cols[factor(round(as.numeric(da$fc),digits=1))]

#copied from https://stackoverflow.com/questions/16875547/using-igraph-how-to-force-curvature-when-arrows-point-in-opposite-directions
autocurve.edges2 <-function (graph, start = 0.5)
{
    cm <- count.multiple(graph)
    mut <-is.mutual(graph)  #are connections mutual?
    el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
        collapse = ":")
    ord <- order(el)
    res <- numeric(length(ord))
    p <- 1
    while (p <= length(res)) {
        m <- cm[ord[p]]
        mut.obs <-mut[ord[p]] #are the connections mutual for this point?
        idx <- p:(p + m - 1)
        if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
            r <- 0
        }
        else {
            r <- seq(-start, start, length = m)
        }
        res[ord[idx]] <- r
        p <- p + m
    }
    res
}

png(paste(fun,"shared_net.png",sep='_'))
par(mfrow=c(2,2))
sapply(1:4,function(x) {
	plot(shared,
	vertex.color=da$cols[da$contrast==unique(da$contrast)[x]],
	vertex.label.cex=1.5,layout=lay,vertex.size=30,
	vertex.label.family="sans",main=names(nets)[x],
	vertex.label.color="black",
	edge.label=c("","-")[factor(as.numeric(ws[,x])<0,
		levels=c("FALSE","TRUE"))],
	edge.label.cex=3,
	edge.width=10*abs(as.numeric(ws[,x])),
	edge.curved=autocurve.edges2(shared))})
dev.off()
#WHAT MAKES SPECIAL THIS EDGES?????????



#MAKES NETWORK TOO LARGE
#↓
#↓
#↓
#########################ADD KNOWN REGULATORY EDGES « net becomes huge
#library(biomaRt)
#files=list.files()
#files=files[grep(fun,files)]
#files=files[grep("moti",files)]
#known=lapply(files,read_tsv,col_names=F)
#known=do.call(rbind,known)
#known=graph.data.frame(known[,1:2])
##focus only on the nodes with shared links 
#known=make_ego_graph(known,
#	nodes=which(V(known)$name%in%V(shared)$name))
#known=as.data.frame(unique(do.call(rbind,lapply(known,get.edgelist))))
#known$type="known"
#final=rbind(cbind(get.edgelist(shared),type="found"),
#			known)
#g=graph.data.frame(final[,1:2])
#
#bps=c("GO:0048568","GO:0061326","GO:0072006","GO:0072073","GO:0098742")
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
#	version=105)
##https://dec2021.archive.ensembl.org
#myannot <- getBM(attributes=c('hgnc_symbol', 
#			'ensembl_gene_id','go_id','name_1006'),
#			filters="ensembl_gene_id",
#			values=V(g)$name,
##filters = 'go', 
#			#values = bps,
#			 mart = mart)
#                   #por gene names y agrupa targets por función?????????
#