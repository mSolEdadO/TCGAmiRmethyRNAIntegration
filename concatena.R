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

i=grep("^rs",rownames(concatenadas$LumA[[2]]),perl=T,invert=T)
#hay que quitar las probes de SNPs
concatenadas$LumA[[2]]=concatenadas$LumA[[2]][i,]
concatenadas$Basal[[2]]=concatenadas$Basal[[2]][i,]
concatenadas$LumB[[2]]=concatenadas$LumB[[2]][i,]
concatenadas$Her2[[2]]=concatenadas$Her2[[2]][i,]
concatenadas$normal[[2]]=concatenadas$normal[[2]][i,]
lapply(concatenadas,function(x) sapply(x,dim))
#$LumA
#     [,1]   [,2]  [,3]
#[1,]  446 395741 13904
#[2,]  331    331   331

#$Basal
#     [,1]   [,2]  [,3]
#[1,]  446 395741 13904
#[2,]  135    135   135

#$LumB
#     [,1]   [,2]  [,3]
#[1,]  446 395741 13904
#[2,]  177    177   177

#$Her2
#     [,1]   [,2]  [,3]
#[1,]  446 395741 13904
#[2,]   75     75    75
#$normal
#     [,1]   [,2]  [,3]
#[1,]  446 395741 13904
#[2,]   75     75    75

save(concatenadas,file="conca/porSubti.RData")
#########################################

