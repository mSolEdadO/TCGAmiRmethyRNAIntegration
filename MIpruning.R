library(igraph)
library(data.table)
###funciones de Memo
readSIF <- function(SIFfile){
   #read a SIF file into R 
   g<-fread(SIFfile, data.table = FALSE)
   g<-g[order(g$V2, decreasing = TRUE),] #ensure ordering
   names(g)[names(g)=="V2"] <- "weight"
   g<-g[,c(1,3,2)]
  return(g)
}
pvalue<-function(mi, n=100){
  alfa = 1.062
  beta = -48.7
  gamma = -0.634
  p = exp(alfa -mi*(-beta + (-gamma * n)))
  return(p)
}

sif=readSIF("parallel-aracne/results/lumAmiRmRNA.sif")
sif$pvalue=pvalue(sif$weight)
sif$q=p.adjust(sif$pvalue)

i=grep("hsa",sif$V1)
j=grep("hsa",sif$V3)
mirs=intersect(i,j)
mrnas=intersect(which(!1:nrow(sif)%in%i),which(!1:nrow(sif)%in%j))
mix=which(!1:nrow(sif)%in%c(mirs,mrnas))
sif$type="mrna-mrna"
sif$type[mirs]="mir-mir"
sif$type[mix]="mrna-mir"

boxplot(q~type,sif)
sapply(c("mir-mir","mrna-mir","mrna-mrna"),function(x) sum(sif$q[sif$type==x]<0.01))
