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
			c("function"=fun,"pair1"=subtys[x],"pair2"=subtys[y],
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
jindx%>%filter(type=="biological process")%>%
ggplot(aes(x=as.numeric(index),y=pair))+geom_violin()+
geom_jitter(height=0.25)+theme(text=element_text(size=16),
	axis.ticks.y=element_blank(),panel.background=element_blank())+
ylab("")+xlab("jaccard index")+geom_vline(xintercept=0.5,
	color="firebrick")
dev.off()
#save it for puma nets
write_tsv(x=jindx,"funcJaccI.tsv")


