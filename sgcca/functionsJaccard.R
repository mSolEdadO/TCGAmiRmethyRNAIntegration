#which features are most enriched?????????????????????
bp=read_tsv("BP.enrichment")
k=read_tsv("KEGG.enrichment")
enriched=list(bp,k)
names(enriched)=c("BP","KEGG")

features=lapply(enriched,function(x) 
	x%>%dplyr::select(subtype,Description,geneID)%>%
	separate_rows(geneID,sep='/',convert=T))

#which functions map the same genes for a pair of datasets?????
i=lapply(features,function(x) 
	x%>%distinct(subtype,Description)%>%count(Description)%>%
	filter(n>1)%>%select(Description)%>%unlist)
#use only functions found in at least 2 datsets
features=lapply(1:2,function(x) 
	features[[x]]%>%filter(Description%in%i[[x]]))

#jaccard index
jacc=function(set1,set2){
	inter=length(intersect(set1,set2))
	un=length(union(set1,set2))
return(inter/un)}
intersect_functions=function(data,fun){
	set=data%>%filter(Description==fun)
	subtys=unique(set$subtype)
	#paired contrast
	mat=do.call(rbind,lapply(1:(length(subtys)-1),function(x) 
		do.call(rbind,lapply((1+x):length(subtys),function(y) 
			c("func"=fun,"pair1"=subtys[x],"pair2"=subtys[y],
				"index"=jacc(set$geneID[set$subtype==subtys[x]],
						   set$geneID[set$subtype==subtys[y]]))))))
return(mat)}
jindx=lapply(features,function(x) 
	data.frame(do.call(rbind,lapply(unique(x$Description),function(y)
		intersect_functions(x,y)))))

#toplot
jindx=do.call(rbind,lapply(1:2,function(x) 
	cbind("type"=c("biological process","KEGG pathway")[x],
		jindx[[x]])))
jindx$pair=paste(jindx$pair1,jindx$pair2)
#jindx%>%count(type)
#                type   n
# biological process 523
#       KEGG pathway  20
#KEGG has to few (and low) points to plot
#plot distributions of Jaccard Index
png("enrichJacc.png")
ggplot(jindx,aes(x=as.numeric(index),y=pair))+
geom_violin()+geom_jitter(height=0.25)+
ggrepel::geom_text_repel(data=jindx[jindx$index>0.5,],
 aes(y=pair,x=index,label=func),hjust=1,nudge_x=1.5,xlim=c(0,1.2),
 direction='y')+theme(text=element_text(size=16),
 axis.ticks.y=element_blank(),panel.background=element_blank(),
 plot.margin=unit(c(.1,5,.1,.1),"cm"))+coord_cartesian(clip="off")+
ylab("")+xlab("jaccard index")+geom_vline(xintercept=0.5,
	color="firebrick")
dev.off()
#save it for puma nets
write_tsv(x=jindx,"funcJaccI.tsv")

##############shared plot like in functionalEnrichment.R
#get functions of interest
i=jindx%>%filter(index>0.5)%>%distinct(func)%>%unlist
#get aaaall realted functions
groups=read_tsv("BP.groups")
#get groups of functions in i
i=groups%>%filter(Description%in%i)%>%distinct(subtype,group)
i=merge(i,groups,by=c("subtype","group"),all.x=T)

#toplot
bp=bp[bp$Description%in%i$Description,]
#count genes & component per function and subtype
edges=bp%>%filter(Description%in%i$Description)%>%
	dplyr::select(subtype,Description,geneID)%>%
	group_by(subtype,Description)
edges=edges%>%group_map(~length(unique(unlist(strsplit(.x$geneID,
	"/")))))%>%unlist%>%cbind(group_keys(edges))
colnames(edges)[1]="genes"
edges1=bp%>%group_by(subtype,Description)%>%tally
colnames(edges1)[3]="components"
edges=merge(edges,edges1,by=c("subtype","Description"),all.x=T)
#add groups info
edges=merge(i,edges,by=c("subtype","Description"),all.y=T)
#trick to repeat figures between datasets
edges$group=edges%>%group_by(subtype)%>%
	group_map(~c(0:7)[factor(.x$group)])%>%unlist
edges$group[is.na(edges$group)]=7
#ADD GSEA INFO
gsea=read_tsv("BP.gsea")
gsea=gsea%>%filter(p.adjust<0.05)

png("BPeshared.png",width=800,height=650)
 ggplot(edges)+geom_point(aes(x=subtype,y=Description,
 	size=components,col=genes,shape=as.character(group)))+
 scale_color_gradient(low="blue",high="red")+
 theme_light(base_size=16)+scale_size(range=c(3,10),trans="log10")+
 theme(axis.ticks=element_blank())+xlab("")+ylab("")+
 scale_shape_manual(values=c(17,18,15,0:2,16),guide="none")+
 annotate("text",
 	y=gsea$Description[gsea$Description%in%edges$Description],
 	x="Basal",label="*",size=7,vjust=.8,hjust=4.2)+
 coord_cartesian(clip="off")
dev.off()
