library(RCy3)
library(igraph)

#set graph object
edges <- read.table("Downloads/4clusters/clust2.tsv",header=T)
g <- graph.data.frame(edges[,1:2],directed=F)
#E(g)$rho <- edges$rho
#E(g)$adjustedPval <- edges$adjustedPval
V(g)$membership <- cluster_louvain(g)$membership
#cytoscape must be running before next line
#g=induced_subgraph(g,1:1000)
gc=createNetworkFromIgraph(g,"myIgraph")
setNodeColorMapping(table.column="membership",
    mapping.type="discrete",
    colors=paletteColorBrewerDark2)

#group by membership
nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
sapply(unique(nodedata$membership),function(x) {
    geneSUIDs <- nodedata$SUID[nodedata$membership==x]
    selectNodes(geneSUIDs, preserve.current.selection = FALSE)
    tag=paste("group",x)
    createGroup(tag)
    collapseGroup(tag)})