library(data.table)
library(ggplot2)
library(gridExtra)
#get data
files=list.files()
sam=files[grep("sort",files)]
ori=files[grep("ori",files)]

observed=lapply(ori,fread)
names(observed)=gsub(".ori",ori)
predi=lapply(sam,fread)
names(predi)=sapply(strsplit(sam,".",fixed=T),
	function(x) x[1])
#join MI values for all the samples per subtype
predi=lapply(names(observed),function(x) predi[names(predi)==x])
predi=lapply(predi,function(x) 
	do.call(cbind,lapply(x,function(y) y$V3)))
#save data
lapply(1:4,function(x) write.table(predi[[x]],
	paste("/home/msoledad/Downloads/",names(observed)[x],
	".sampled",sep=''),sep='\t',quote=F,
	row.names=paste(observed[[x]]$V1,observed[[x]]$V3),col.names=F))
#get stuff for zeta score
statistics=lapply(predi,function(x) 
	cbind(mean=apply(x,1,mean),
          sd=apply(x,1,sd)))
#compare sampled vs complete
paplot=lapply(1:4,function(x) 
	data.frame(cbind(sampled=statistics[[x]][,1],
		complete=observed[[x]]$V2,
		enriched=observed[[x]]$enriched),stringsAsFactors=F))
paplot=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(names(observed)[x],paplot[[x]]))))
colnames(paplot)[1]="subtype"
png("/home/msoledad/Downloads/stability.png",width=800)
ggplot(paplot,aes(x=sampled,y=complete,color=factor(enriched)))+
geom_point()+facet_wrap(~subtype)+theme(text=element_text(size=18))+
scale_color_discrete(name="biological process",
	labels=c("FALSE","TRUE"))
dev.off()
#get z-score
zetas=lapply(1:4,function(x)
	(observed[[x]]$V2-statistics[[x]][,1])/statistics[[x]][,2])
zetas=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(names(observed)[x],zetas[[x]]))))
colnames(zetas)=c("subtype","z")
zetas$z=as.numeric(as.character(zetas$z))
#plot z distri
p1=ggplot(zetas,aes(x=z))+geom_density(aes(fill=subtype,
	color=subtype,y=..scaled..),alpha=0.3)+xlim(-1,10)+
	theme(text=element_text(size=18),legend.position="n")+
	scale_color_manual(values=scales::hue_pal()(5)[c(1,3:5)])+
	scale_fill_manual(values=scales::hue_pal()(5)[c(1,3:5)])
p2=ggplot(zetas,aes(x=subtype,y=z,color=subtype))+geom_boxplot()+
	theme(text=element_text(size=18))+xlab("")+
	scale_color_manual(values=scales::hue_pal()(5)[c(1,3:5)])
png("/home/msoledad/Downloads/zetas.png",width=900)
grid.arrange(p1,p2,ncol=2)
dev.off()
