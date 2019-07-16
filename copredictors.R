
regus=lapply(sifs2,function(x) sapply(unique(x[,4]),function(y) x[x[,4]==y,5]))
Jindex=lapply(regus,function(z){
	distancias=matrix(ncol=length(z),nrow=length(z),dimnames=list(rownames=names(z),colnames=names(z)))
	i=which(upper.tri(distancias),arr.ind=T)
	i=i[order(i[,1]),]
	distancias[i]=unlist(lapply(1:(length(z)-1),function(x) sapply((x+1):length(z),function(y) length(intersect(z[[x]],z[[y]]))/length(c(z[[x]],z[[y]])))))
	i=which(lower.tri(distancias),arr.ind=T)
	i=i[order(i[,1]),]
	distancias[i]=distancias[upper.tri(distancias)]
	diag(distancias)=1
	return(distancias)})
edgelists=lapply(Jindex,function(x) cbind(rownames(x)[which(x>0,arr.ind=T)[,1]],colnames(x)[which(x>0,arr.ind=T)[,2]],x[which(x>0)]))
g=lapply(edgelists,function(x) {g=graph.edgelist(x[,1:2],directed=F);E(g)$weight=1-as.numeric(x[,3]);as.matrix(g[unique(x[,1]),unique(x[,2])])})
g$Basal[g$Basal==0]=NA
g$Her2[g$Her2==0]=NA
g$LumA[g$LumA==0]=NA
g$LumB[g$LumB==0]=NA
g$normal[g$normal==0]=NA

pdf("copredictors.pdf")
> heatmap.2(g$Basal,col=rev(heat.colors(100)),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$Basal)]],colCol=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%colnames(g1$Basal)]],key=F,main="Basal",na.color="black",Rowv=NA,Colv=NA,dendrogram='n',trace="b",margins=c(4,6))
> heatmap.2(g$Her2,col=rev(heat.colors(100)),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$Her2)]],colCol=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%colnames(g1$Her2)]],key=F,main="Her2",na.color="black",Rowv=NA,Colv=NA,dendrogram='n',trace="b",margins=c(4,6))
> heatmap.2(g$LumA,col=rev(heat.colors(100)),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$LumA)]],colCol=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%colnames(g1$LumA)]],key=F,main="LumA",na.color="black",Rowv=NA,Colv=NA,dendrogram='n',trace="b",margins=c(4,6))
> heatmap.2(g$LumB,col=rev(heat.colors(100)),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$LumB)]],colCol=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%colnames(g1$LumB)]],key=F,main="LumB",na.color="black",Rowv=NA,Colv=NA,dendrogram='n',trace="b",margins=c(4,6))
> heatmap.2(g$normal,col=rev(heat.colors(100)),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$normal)]],colCol=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%colnames(g1$normal)]],key=F,main="normal",na.color="black",Rowv=NA,Colv=NA,dendrogram='n',trace="b",margins=c(4,6))
> dev.off()
#########################################################3
g=lapply(sifs2,function(x) {g=graph.edgelist(x[,4:5],directed=F);E(g)$weight=as.numer
ic(x[,3]);as.matrix(g[unique(x[,4]),unique(x[,5])])})
g=lapply(g,function(x) x[,colSums(x>0)>1])
g=lapply(g,function(x) x[rowSums(x==0)<ncol(x),])
g1=g
g1$Basal[g1$Basal==0]=NA
g1$Her2[g1$Her2==0]=NA
g1$LumA[g1$LumA==0]=NA
g1$LumB[g1$LumB==0]=NA
g1$normal[g1$normal==0]=NA
pdf("copredictors.pdf")
heatmap.2(g1$Basal,scale='n',col=colorRampPalette(colors = c("red4","red"))(10),na.color="black",dendrogram='n',Colv=NA,Rowv=NA,symm=F,breaks=seq(0.01,0.1,length=11),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$Basal)]],margins=c(8,5),key=F,trace='n')


heatmap.2(g1$Her2,scale='n',col=colorRampPalette(colors = c("red4","red"))(10),na.color="black",dendrogram='n',Colv=NA,Rowv=NA,symm=F,breaks=seq(0.01,0.1,length=11),colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$Her2)]],margins=c(10,10),key=F,trace='n')

heatmap.2(g1$LumA,scale='n',col=colorRampPalette(colors = c("green","green4","red4","red"))(10),na.color="black",dendrogram='n',Colv=NA,Rowv=NA,symm=F,colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$LumA)]],margins=c(8,5),trace='n',breaks=c(seq(-0.5,-0.1,length=5),0,seq(0.01,0.5,length=5)),key=F)

heatmap.2(g1$LumB,scale='n',col=colorRampPalette(colors = c("green","green4","red4","red"))(10),na.color="black",dendrogram='n',Colv=NA,Rowv=NA,symm=F,colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$LumB)]],margins=c(8,5),trace='n',breaks=c(seq(-0.5,-0.1,length=5),0,seq(0.01,0.5,length=5)),key=F)

heatmap.2(g1$normal,scale='n',col=colorRampPalette(colors = c("green","green4","red4","red"))(10),na.color="black",dendrogram='n',Colv=NA,Rowv=NA,symm=F,colRow=ggsci::pal_jama()(4)[pam50$class[pam50$hgnc_symbol%in%rownames(g1$normal)]],margins=c(8,5),trace='n',breaks=c(seq(-0.5,-0.1,length=5),0,seq(0.01,0.5,length=5)),key=F)
> dev.off()
##############################
i=names(which(table(unlist(lapply(g,colnames)))>1))
g=sapply(g,function(x) dim(x[,colnames(x)%in%i]))
total=as.matrix(rbind.fill(lapply(g,as.data.frame)))
rownames(total)=unlist(sapply(g,rownames))


