load("ini/subtiTMMArsyn.RData")
load("ini/mirTMMARSyn.RData")
load("ini/imptMethy.RData")

shared=lapply(c("LumA","Basal","LumB","Her2","normal"),function(x) sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
	substr(colnames(y[[which(names(y)==x)]]),1,12)))
shared=sapply(1:5,function(x) intersect(intersect(shared[[x]][[1]],shared[[x]][[2]]),shared[[x]][[3]]))
names(shared)=c("LumA","Basal","LumB","Her2","normal")
concatenadas=lapply(1:5,function(x) sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
	y[[which(names(y)==names(shared)[x])]][,substr(colnames(y[[which(names(y)==names(shared)[x])]]),1,12)%in%shared[[x]]]))
names(concatenadas)=names(shared)
concatenadas=lapply(concatenadas,function(x) sapply(x,function(y) y[,order(substr(colnames(y),1,12))]))
save(concatenadas,normal,file="conca/porSubti.RData")
#########################################
library(pbapply)
library(mixOmics)

escaladas=pblapply(concatenadas,function(x) lapply(x,function(y) scale(y,center=T,scale=T)))
vari=pblapply(escaladas,function(x) sapply(x,var))
pdf("varConcate.pdf")
par(mfrow=c(4,3))
lapply(1:4,function(x) sapply(1:3,function(y) plot(density(vari[[x]][,y]),main=paste(names(concatenadas)[x],c("miRNA","methy","mRNA")[y],sep=':'))))
dev.off()
#vari1=pblapply(escaladas,function(x) sapply(t(x),var))

#indi.pca=lapply(concatenadas,function(x) pca(t(x),ncomp=50,center=T,scale=T))
pdf("varCum.pdf")
par(mfrow=c(4,3))
lapply(1:5,function(x) lapply(1:3,function(y) 
	plot(indi.pca[[x]][,y]$cum.var,main=paste(names(concatenadas)[x],c("miRNA","methy","mRNA")[y],sep=':'),ylim=c(0,1),ylab="cum.var")
	abline(h=50,col="red")))
dev.off()
#lapply(indi.pca,function(x) apply(x,2,function(y) tail(y$cum.var,3)))
#$LumA
#        miRNA      methy      mRNA
#PC98  0.6358782 0.5697302 0.5798319
#PC99  0.6390467 0.5728921 0.5832214
#PC100 0.6421894 0.5760480 0.5865936
#$Basal
#           [,1]      [,2]      [,3]
#PC98  0.8962268 0.8873851 0.8678912
#PC99  0.9001018 0.8915093 0.8725604
#PC100 0.9038882 0.8956038 0.8771721
#$LumB
#           [,1]      [,2]      [,3]
#PC98  0.7988851 0.7451324 0.7511437
#PC99  0.8027174 0.7495495 0.7555117
#PC100 0.8065252 0.7539502 0.7598363
#$Her2
#           [,1]      [,2]      [,3]
#PC98  0.9336252 0.9537041 0.9153974
#PC99  0.9375400 0.9565811 0.9203258
#PC100 0.9414300 0.9593893 0.9252081

concatenadas1=lapply(concatenadas,function(x) do.call(rbind,x))
probes=as.data.frame(cbind(rownames(concatenadas1$LumA),c(rep("miRNA",1588),rep("methy",395808),rep("mRNA",13904))))
colnames(probes)=c("name","type")
concate.pcas=pblapply(concatenadas1,function(x) pca(t(x),ncomp=10,center=T,scale=T))
#sapply(1:4,function(x) which(concate.pcas[[x]][[15]]>0.5)[1])
#<NA> PC37 <NA> PC29 
#  NA   37   NA   29 
png("subti10compo.png",width=200,height=800)
par(mfrow=c(4,1))
sapply(1:4,function(x) plot(concate.pca[[x]]$cum.var,main=names(concatenadas)[x],ylim=c(0,1),ylab="cum.var")
	dev.off()
#convert +append pcaIndi.png subti10compo.png lala.png

png("Her2dataCorr.png")
plotVar(concate.pcas[[4]],comp=c(1,2),var.names=F,title="Her2",legend=T,
	col=list(rep(c("cornflowerblue","darkmagenta","brown1"),table(probes$type[probes$name%in%rownames(concatenadas1[[4]])]))))
#legend("top",legend=c("methy","miRNA","mRNA"),fill=c("cornflowerblue","green","brown1"))
dev.off()
png("LumAdataCorr.png")
plotVar(concate.pcas[[1]],comp=c(1,2),var.names=F,title="LumA",legend=T,
	col=list(rep(c("cornflowerblue","darkmagenta","brown1"),table(probes$type))))
dev.off()
png("LumBdataCorr.png")
plotVar(concate.pcas[[2]],comp=c(1,2),var.names=F,title="LumB",legend=T,
	col=list(rep(c("cornflowerblue","darkmagenta","brown1"),table(probes$type))))
dev.off()
png("BasaldataCorr.png")
plotVar(concate.pcas[[3]],comp=c(1,2),var.names=F,title="Basal",legend=T,
	col=list(rep(c("cornflowerblue","darkmagenta","brown1"),table(probes$type))))
dev.off()