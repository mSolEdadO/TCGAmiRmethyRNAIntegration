library(data.table)
library(mixOmics)

raw=fread("Her2.mtrx")
raw=as.matrix(raw[,2:ncol(raw)],rownames=raw$V1)
raw=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,function(x) raw[x[1]:x[2],])
normi=fread("Her2.normalized")
normi=as.matrix(normi[,2:ncol(normi)],rownames=normi$V1)
normi=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,function(x) normi[x[1]:x[2],])
names(raw)=c("CpGs","transcripts","miRNAs")
names(normi)=c("CpGs","transcripts","miRNAs")

#note: the penalty parameters will need to be tuned
#op1:cca=wrapper.sgcca(raw,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op2:cca.n=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=F)
#op3:cca.n1=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op1 y op3 seleccionan las mismas features,
#pero con algunos loadgings y loadings.star distintos 
#op1 converge (cca$crit) un poco antes que op3
# NO tiene caso normalizar por eigenvalue
#selectVar() se basa en loadings != 0, NO en loadings.star

#checa tune.block.splsda para implementar con caret
#y tune.spls con X concatenada y escalada?
