#!/usr/bin/env Rscript 

#befores
#concatenadas=lapply(concantenas,function(x) do.call(rbind,x))
#lapply(1:5,function(x) write.table(concatenadas[[x]],
#	paste("parallel-aracne/data/",names(concatenadas)[x],".txt",
#		sep=''),sep='\t',quote=F))

library("data.table")
library(pbapply)
  
#subtype = commandArgs(trailingOnly=TRUE)
subtype="normal"
#print("Cargo datos")
data=fread(paste("parallel-aracne/data/",subtype,".txt",sep=""))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#write.table(top$normal,"normal.temp",sep="\t",quote=F,row.names=F,col.names=c("node1","MI","node2","enriched"))
MItop=fread(paste("parallel-aracne/",subtype,".temp",sep=""))
#print("Calculos")
#get correlation and p-value for given pairs
corres=pbapply(MItop,1,function(x) 
	unlist(cor.test(data[rownames(data)==x[1]],
		data[rownames(data)==x[3]])[3:4]))
MItop=cbind(MItop,t(corres))
MItop$q=p.adjust(MItop$p.value,"fdr")
colnames(MItop)[5:6]=c("pval","cor")

write.table(MItop,paste("parallel-aracne/",subtype,".temp",sep=""),
	sep='\t',quote=F,row.names=F)

################bring all together
files=list.files("Downloads",full.names=T)
files=files[grep("temp",files)]
corres=lapply(files,fread)
names(corres)=gsub("Downloads/","",gsub(".temp","",files))

plots=lapply(1:5,function(x) {
 ggplot(corres[[x]],aes(x=MI,y=cor))+geom_point(aes(color=enriched))+
 geom_ribbon(aes(ymin=min(cor[!pval<0.05]),
 	ymax=max(cor[!pval<0.05])),alpha=0.3)+
 scale_x_continuous(expan=c(0,0))+
 theme(text=element_text(size=18),legend.position="bottom")+
 scale_colour_manual(name="biological process",
 	values=c("#00BFC4","#F8766D"))})

png("Downloads/Basal_correlation.png")
 plots[[1]]
dev.off()
png("Downloads/Her2_correlation.png")
 plots[[2]]
dev.off()
png("Downloads/LumA_correlation.png")
 plots[[3]]
dev.off()
png("Downloads/LumB_correlation.png")
 plots[[4]]
dev.off()
png("Downloads/normal_correlation.png")
 plots[[5]]
dev.off()

#how many from each kind
temp=lapply(corres,function(x) x[x$pval>0.05,c(1,3,4)])
lala=lapply(temp,function(x) 
	table(data.frame(cbind(paste(substr(x$node1,1,1),
		substr(x$node2,1,1)),x$enriched))))
lolo=lapply(corres,function(x) 
	table(data.frame(cbind(paste(substr(x$node1,1,1),
		substr(x$node2,1,1)),x$enriched))))