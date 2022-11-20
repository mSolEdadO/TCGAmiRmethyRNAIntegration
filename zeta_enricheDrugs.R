library(tidyverse)
library(ggplot2)

real=read_tsv("BreastenricheDrugs.tsv")
real$pathway=gsub("comm",":comm",
	gsub("clust","class",real$pathway))

toplot=real%>%select(pert_iname,pathway,NES)%>%
group_by(pathway,pert_iname)%>%summarise(max=max(NES),min=min(NES))
#plot positive NES
png("NES+.png")
toplot%>%filter(max>0)%>%ungroup%>%complete(pathway,pert_iname)%>%
ggplot(aes(pathway,pert_iname))+geom_tile(aes(fill=max))+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
	axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_fill_viridis_c(option="magma")+
xlab("MoNet Cluster:Community")+ylab("Drugs")+labs(fill="NES")+
ggtitle("Normalized Enrichment Score
	(maximum positive value across experiments)")
dev.off()
#plot negative
png("NES-.png")
toplot%>%filter(min<0)%>%ungroup%>%complete(pathway,pert_iname)%>%
ggplot(aes(pathway,pert_iname))+geom_tile(aes(fill=min))+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
	axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_fill_viridis_c(option="magma",direction=-1)+
xlab("MoNet Cluster:Community")+ylab("Drugs")+labs(fill="NES")+
ggtitle("Normalized Enrichment Score
	(minimum negative value across experiments)")
dev.off()

files=list.files("92742",full.names=T)
files=c(files,list.files("70138",full.names=T))
random=lapply(files,read_tsv)
random=do.call(rbind,lapply(1:length(files),function(x)
 cbind(gsub(".gsea","",sapply(strsplit(files,"/"),function(y) y[2]))[x],
 	random[[x]])))
#[1] 17108000
random=random[random$sig_id%in%real$sig_id,]
#[1] 962000
zscore=function(real_id,real_pathway,pval){
	pvalues=random%>%
	filter(sig_id==real_id&pathway==real_pathway)%>%
	select(padj)%>%unlist
	pvalues=pvalues[!is.na(pvalues)]
	m=mean(c(pval,pvalues))
	d=sd(c(pval,pvalues))
	z=(pvalues-m)/d
	return(c(zscore=(pval-m)/d,min=min(z),max=max(z)))
}
zetas=apply(real,1,function(x) 
	zscore(x["sig_id"],x["pathway"],as.numeric(x["padj"])))
real=cbind(real,t(zetas))
write_tsv(real,"BreastenricheDrugs-z.tsv")

toplot=real%>%select(pert_iname,pathway,zscore)%>%
group_by(pathway,pert_iname)%>%
summarise(meanz=mean(zscore),maxz=max(zscore))
#toplot=toplot %>% ungroup%>%complete(pathway, pert_iname)
png("drug_meanZscore.png")
toplot%>%filter(meanz<(-3))%>%ungroup%>%
complete(pathway,pert_iname)%>%
ggplot(aes(pathway,pert_iname))+geom_tile(aes(fill=meanz))+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
	axis.ticks.x=element_blank(),
	axis.ticks.y=element_blank())+
scale_fill_viridis_c(option="magma",direction=-1)+
xlab("MoNet Cluster:Community")+ylab("Drugs")+labs(fill="zscore")+
ggtitle("Mean zscore across experiments
	Threshold = -3")
dev.off()
png("drug_maxZscore.png")
toplot%>%filter(maxz<(-3))%>%ungroup%>%
complete(pathway,pert_iname)%>%
ggplot(aes(pathway,pert_iname))+
geom_tile(aes(fill=maxz))+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
	axis.ticks.x=element_blank(),
	axis.ticks.y=element_blank())+
scale_fill_viridis_c(option="magma",direction=-1)+
xlab("MoNet Cluster:Community")+ylab("Drugs")+labs(fill="zscore")+
ggtitle("Maximum zscore across experiments
	Threshold = -3")
dev.off()
