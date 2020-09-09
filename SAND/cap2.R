library(igraph)
library(igraphdata)

data(Ecoli.data)
g=graph.adjacency(regDB.adj)
#matrix to ~sif
get.edgelist(g)
#########tamaño de la red
length(V(g))
#[1] 153
vcount(g)#hace lo mismo
#[1] 153
length(E(g))
#[1] 220
ecount(g)#hace lo mismo
#[1] 220
summary(g) #da lo mismo
#IGRAPH 9f4fb1e DN-- 153 220 -- 
print_all(g)
#print all edges?219#################################
str(g)
#qué hace?

is.weighted(g)
#[1] FALSE
is.directed(g)
#[1] TRUE

#muchos nodos tienen grado 0
sum(degree(g, mode="out")>0)
#[1] 100
sum(degree(g, mode="in")>0)
#[1] 63
#tiene varios componentes
is.connected(g)
#[1] FALSE
is.connected(g, mode="weak")
#[1] FALSE
clusters(g)$no
#[1] 39

diameter(g)
#[1] 7 es poco? qué quiere decir?

V(g)$color <- "red"
plot(g)
list.vertex.attributes(g)
#[1] "name"  "color"
#################subredes
#me quedo con 50 nodos de la original
h1=induced.subgraph(g,sample(1:V(g),50))
summary(h1)
#IGRAPH 470dcbe DN-- 50 14 -- 

#quito 50 nodos de la original
h2 <- g - vertices(sample(1:length(V(g)),50))
summary(h2)
#IGRAPH 1de63c7 DN-- 103 101 -- 

#le agrego 4 aristas a la de 50 nodos
h1.1 <- h1 + edges(c(4,6),c(4,7),c(5,6),c(6,7))
summary(h1.1)
#IGRAPH 6ccba78 DN-- 50 18 -- 

#le agrego 10 nodos
h1.1 <- h1.1 + vertices(sample(1:length(V(g)),10))
summary(h1.1)
#IGRAPH ab4d9d3 DN-- 60 18 -- 

summary(graph.union(h1,h2))
#IGRAPH 8baeb8b DN-- 123 109 -- 
summary(graph.difference(h1,h2))
#IGRAPH 2e47f9b DN-- 50 8 --
summary(graph.intersection(h1,h2))
#IGRAPH 75d7578 DN-- 123 6 -- ################???????????????????????

##########multi-graph?
is.simple(g)
#[1] TRUE
mg <- g + edge(1,1)
is.simple(mg)
#[1] FALSE
g1=simplify(mg)
is.simple(g1)
#[1] TRUE
############neighbors
sapply(1:vcount(h1),function(x) neighbors(g,x))
############grado
