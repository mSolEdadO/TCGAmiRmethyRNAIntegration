#msoledad@castillo:~/parallel-aracne$ ;while [ ! -f lumAmiRmRNA.14350.hsa-mir-99b.adj ]; do condor_release -all;sleep 120; done

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
pvalue<-function(mi, n){
  alfa = 1.062
  beta = -48.7
  gamma = -0.634
  p = exp(alfa -mi*(-beta + (-gamma * n)))
  return(p)
}

sif=readSIF("parallel-aracne/results/lumAmiRmRNA.sif")

i=grep("hsa",sif$V1)
j=grep("hsa",sif$V3)
mirs=intersect(i,j)
mrnas=intersect(which(!1:nrow(sif)%in%i),which(!1:nrow(sif)%in%j))
mix=which(!1:nrow(sif)%in%c(mirs,mrnas))
sif$type="mrna-mrna"
sif$type[mirs]="mir-mir"
sif$type[mix]="mrna-mir"
table(sif$type)
#  mir-mir  mrna-mir mrna-mrna 
#    99235   6201184  96653656 
png("mi1.png")
ggplot(sif,aes(x=weight))+geom_density(aes(group=type,colour=type,fill=type),alpha=0.3)+xlim(0.25,1)
dev.off()


sif$pvalue=pvalue(sif$weight,331)
sif$q=p.adjust(sif$pvalue,"fdr")
sapply(unique(sif$type),function(x) sum(sif$q[sif$type==x]<0.01))
#  mir-mir mrna-mrna  mrna-mir 
#     1976     29972        18 
boxplot(q~type,sif)
