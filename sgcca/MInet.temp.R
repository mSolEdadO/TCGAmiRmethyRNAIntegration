library(igraph)
library(tidyverse)
library(biomaRt)
library(RCy3)


edges=read.table("LumB.sif")
g=graph.data.frame(edges,directed=F)
E(g)$MI=edges$V3
g=igraph::simplify(g,edge.attr.comb="first")

V(g)$type=substr(V(g)$name,1,1)
tfs=read_tsv("../../TFtargets.tsv")
tfs=unique(unlist(strsplit(tfs$TF,',')))
V(g)$type[V(g)$name%in%tfs]="TF"
targets=readLines("node.list.temp")
V(g)$name[V(g)$type=="E"][!V(g)$name[V(g)$type=="E"]%in%targets]
#[1] "ENSG00000115590" "ENSG00000082701"
g=delete_vertices(g,c("ENSG00000115590","ENSG00000082701"))
targets=c(targets,unlist(sapply(targets,function(x) neighbors(g,which(V(g)$name==x))$name)))
g=delete_vertices(g,V(g)$name[!V(g)$name%in%targets])

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
"entrezgene_accession"), mart=mart)
V(g)$symbol=V(g)$name
myannot=myannot[myannot$ensembl_gene_id%in%V(g)$name,]
V(g)$symbol[V(g)$name%in%myannot$ensembl_gene_id]=myannot$entrezgene_accession[order(match(myannot$ensembl_gene_id,V(g)$name))]

methy=read_tsv("../../MapMethy.tsv")
methy=methy[methy$IlmnID%in%V(g)$name,]
methy$UCSC_RefGene_Name=sapply(strsplit(methy$UCSC_RefGene_Name,";"),function(x) paste(unique(x),collapse=";"))
methy$UCSC_RefGene_Group=sapply(strsplit(methy$UCSC_RefGene_Group,";"),function(x) paste(unique(x),collapse=";"))
methygroups=methy%>%group_by(UCSC_RefGene_Name,UCSC_RefGene_Group)%>%group_map(~unique(.x$IlmnID))

V(g)$symbol[V(g)$name%in%methy$IlmnID]=methy$UCSC_RefGene_Name[order(match(methy$IlmnID,V(g)$name))]
V(g)$symbol[V(g)$symbol=="NA"]=V(g)$name[V(g)$symbol=="NA"]
#write_graph(g,"LumB-endocrineResist.graphml",format="graphml")
V(g)$size=degree(g)

gc=createNetworkFromIgraph(g,"myIgraph")
setNodeShapeDefault('ELLIPSE')
setNodeColorMapping(table.column="type",
    mapping.type="discrete",
    colors=RColorBrewer::brewer.pal(n = 4, name = "Dark2"))
setNodeSizeMapping(table.column="size")
setNodeLabelMapping(table.column="symbol")
setEdgeLineWidthMapping(table.column="MI")

nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
i=nodedata$symbol[nodedata$type=='c']
i=unique(i[duplicated(i)])
sapply(i,function(x) {
    geneSUIDs <- nodedata$SUID[nodedata$symbol==x&nodedata$type=="c"]
    selectNodes(geneSUIDs, preserve.current.selection = FALSE)
    tag=paste("group",x)
    createGroup(tag)
    collapseGroup(tag)})

edges$type=substr(edges$V2,1,1)
library(ggplot2)
ggplot(edges,aes(x=V3))+geom_density(aes(fill=type),alpha=0.3)
