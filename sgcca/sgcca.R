#plots=lapply(levels(results$omic),function(i) 
#	ggplot(results[results$omic==i,],aes(x=nfeatures,
#		y=AVE,color=sparsity))+geom_point()+ggtitle(i)+
#		scale_x_continuous(trans="log10")+
#		theme(text=element_text(size=18)))
#png(paste(subtype,BP,"png",sep='.'))
#	print({grid.arrange(plots[[1]],plots[[2]],plots[[3]])})
#no facet_wrap coz u want indy axes
#dev.off()

#turn it analyticalish by choosing the largest fall in AVE
sc=which.max(sapply(14:2,function(x) (results$AVE[x]-results$AVE[x-1])))
sc=results$sparsity[14:2][sc]
st=which.max(sapply(28:16,function(x) (results$AVE[x]-results$AVE[x-1])))
st=results$sparsity[28:16][st]
sm=which.max(sapply(42:30,function(x) (results$AVE[x]-results$AVE[x-1])))
sm=results$sparsity[42:30][sm]

##############################GET SELECTED FEATURES
library(igraph)

print("Final SGCCA")
temp=wrapper.sgcca(X=data,penalty=c(sc,sm,st,1),scale=T,
	scheme="centroid",ncomp=components)#ncomp to explain 50% of transcripts matrix according to mfa.R
output=rbind(rowSums(do.call(rbind,temp$AVE$AVE_X)),temp$penalty)
rownames(output)=c("sum_AVE","penalty")
write.table(output,paste(subtype,BP,"descriptors",sep='.'),sep='\t',
	quote=F)
#selectVar() is the same than loadings != 0, NOT loadings.star != 0
#output=lapply(temp$loadings,function(x) x[rowSums(x!=0)>0,])
#lapply(1:4,function(x) 
#	write.table(output[[x]],
#		paste(subtype,BP,names(output)[x],sep='.'),
#		sep='\t',quote=F))

source("function_networkAlt.R")#with no plotting
g=lapply(1:components,function(x) network(temp,comp=list(CpGs=x,
	TFs=x,miRNAs=x,BP=x),blocks=1:4)$gR)#to avoid pdf issues
edges=do.call(rbind,lapply(g,function(x) 
	as.data.frame(cbind(get.edgelist(x),E(x)$weight))))
colnames(edges)=c("source","target","corr")
write.table(edges,paste(subtype,BP,"net",sep='.'),sep='\t',quote=F,
	row.names=F)
#####ISSUE
#como tienes un número diferente de features para cada BP, no sabes
#si un BP tiene mejor AVE que otros, por los valores de penalty
#o porque de hecho el BP no está tan asociado con las omicas
#si fijas los penalties, igualmente las diferencias del AVE pueden
#deberse a que los penalities no son apropiados para el BP

#puedes ¿debes? repetir con los penalties ajustados y 
#matrices aleatorizadas 