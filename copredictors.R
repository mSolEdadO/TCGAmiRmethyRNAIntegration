
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
