library(data.table)
library(igraph)

mi=fread("GO0000070.sif",skip=2)
cor=fread("descri/Her2.GO:0000070.net")
mi=mi[,c(1,3,2)]
g.mi=graph.data.frame(mi[,1:2],directed=F)
E(g.mi)$MI=as.numeric(mi$V2)
#clear repeated edges
g.mi=simplify(g.mi,edge.attr.comb="first")
mi=as.data.frame(get.edgelist(g.mi))
mi$MI=E(g.mi)$MI


pnlty.change=wrapper.sgcca(X=data,penalty=c(0.01,0.05,0.1,1),scale=T,scheme="centroid",ncomp=15)
g=lapply(1:15,function(x) network(pnlty.change,comp=list(CpGs=x,TFs=x,miRNAs=x,BP=x),blocks=1:4)$gR)
pnlty.change=do.call(rbind,lapply(g,function(x) 
as.data.frame(cbind(get.edgelist(x),E(x)$weight))))
pnlty.g=graph.data.frame(edges[,1:2],directed=F)

desi
       CpGs miRNAs TFs BP
CpGs      0      0   0  1
miRNAs    0      0   0  1
TFs       0      0   0  1
BP        1      1   1  0
desi.change=wrapper.sgcca(X=data,penalty=c(0.02,0.06,0.2,1),scale=T,scheme="centroid",ncomp=15,design=desi)
g=lapply(1:15,function(x) network(desi.change,comp=list(CpGs=x,TFs=x,miRNAs=x,BP=x),blocks=1:4)$gR)
desi.change=do.call(rbind,lapply(g,function(x) as.data.frame(cbind(get.edgelist(x),E(x)$weight))))
desi.g=graph.data.frame(desi.change[,1:2],directed=F)
sapply(list(ori.g,pnlty.g,desi.g),vcount)
[1] 4627 1202 4735
sapply(list(ori.g,pnlty.g,desi.g),ecount)
[1] 121506  14285 114188
intersection(ori.g,desi.g,keep.all.vertices=F)
IGRAPH ce6c724 UN-- 866 8940 -- 
intersection(ori.g,pnlty.g,keep.all.vertices=F)
IGRAPH 699f806 UN-- 765 7962 -- 
