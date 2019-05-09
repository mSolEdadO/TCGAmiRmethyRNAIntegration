#build a sif from adjs by chunk
adjs=list.files()
adjs=adjs[3:36]
adjs=lapply(adjs,function(x) readLines(x))
adjs=lapply(adjs,function(y) sapply(y,function(x) unlist(strsplit(x,"\t"),use.names=F)))
adjs=as.data.frame(do.call(rbind,sapply(adjs,function(y) 
  do.call(rbind,sapply(y,function(x) 
    cbind(x[1],t(matrix(as.character(x[2:length(x)]),nrow=2))))))))
colnames(adjs)=c("s","t","mi")
i=paste(adjs$s,adjs$t,sep='-')
dupis=which(duplicated(i))
#adjs[dupis,]
adjs=adjs[!dupis,]

#classify links
mirnas=list(grep("hsa",adjs$s),grep("hsa",adjs$t))
mrnas=list(grep("ENSG",adjs$s),grep("ENSG",adjs$t))
cpgs=list((1:nrow(adjs))[!(1:nrow(adjs))%in%c(mirnas[[1]],mrnas[[1]])],(1:nrow(adjs))[!(1:nrow(adjs))%in%c(mirnas[[2]],mrnas[[2]])])
adjs$type="cpg-mrna"
adjs$type[intersect(mirnas[[1]],mirnas[[2]])]="mir-mir"
adjs$type[intersect(mrnas[[1]],mrnas[[2]])]="mrna-mrna"
adjs$type[intersect(cpgs[[1]],cpgs[[2]])]="cpg-cpg"
adjs$type[c(intersect(mirnas[[1]],mrnas[[2]]),intersect(mrnas[[1]],mirnas[[2]]))]="mir-mrna"
adjs$type[c(intersect(mirnas[[1]],cpgs[[2]]),intersect(cpgs[[1]],mirnas[[2]]))]="mir-cpg"

#describe nodes by type
table(adjs$type)
#  cpg-cpg  cpg-mrna   mir-cpg   mir-mir  mir-mrna mrna-mrna 
# 29418120      7811     44248       526        22      9735 
sapply(unique(adjs$type),function(x) summary(as.numeric(as.vector(adjs$mi[adjs$type==x]))))
#        cpg-cpg mir-cpg cpg-mrna mrna-mrna mir-mir mir-mrna
#Min.     0.1797  0.1797   0.1797    0.1797  0.1798   0.1806
#1st Qu.  0.1936  0.1902   0.1834    0.1914  0.2029   0.1870
#Median   0.2137  0.2055   0.1890    0.2086  0.2421   0.1919
#Mean     0.2295  0.2146   0.1948    0.2229  0.2996   0.1980
#3rd Qu.  0.2488  0.2307   0.1994    0.2393  0.3465   0.2003
#Max.     0.6574  0.4227   0.4173    0.5552  0.6700   0.2484
png("mi1.png")
ggplot(adjs,aes(x=mi))+geom_density(aes(group=type,colour=type,fill=type),alpha=0.3)
dev.off()
