load("ini/subtiTMMArsyn.RData")
load("ini/mirTMMARSyn.RData")
load("ini/imptMethy.RData")

shared=lapply(c("LumA","Basal","LumB","Her2","normal"),function(x) 
			sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
				substr(colnames(y[[which(names(y)==x)]]),1,12)))
shared=sapply(1:5,function(x) 
			intersect(intersect(shared[[x]][[1]],shared[[x]][[2]]),shared[[x]][[3]]))
names(shared)=c("LumA","Basal","LumB","Her2","normal")
concatenadas=lapply(1:5,function(x) sapply(list(mirSubti,Msubti,TMMArsyn),function(y) 
	y[[which(names(y)==names(shared)[x])]][,substr(colnames(y[[which(names(y)==names(shared)[x])]]),1,12)%in%shared[[x]]]))
names(concatenadas)=names(shared)
concatenadas=lapply(concatenadas,function(x) sapply(x,function(y) y[,order(substr(colnames(y),1,12))]))
sapply(concatenadas,function(x) ncol(x[[1]]))
#  LumA  Basal   LumB   Her2 normal 
#   331    135    177     75     75 

save(concatenadas,methyDesign,subtipos,file="porSubti.RData")
#########################################

