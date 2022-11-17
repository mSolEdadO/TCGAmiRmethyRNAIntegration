library(tidyverse)
library(igraph)
library(RCy3)

bp=read_tsv("BP.enrichment")
k=read_tsv("KEGG.enrichment")
enriched=list(BP=bp,KEGG=k)

##############FIND "CROSSLINKED" FUNCTIONS
coenriched=do.call(rbind,enriched)%>%group_by(subtype)%>%
	group_map(~table(.x[,c("Description","component")]))
#matrix with the component intersections
intersection=lapply(coenriched,function(x) crossprod(t(x)))
#matrix with the component unions
union=lapply(coenriched,function(z) sapply(1:nrow(z),function(x) 
	sapply(1:nrow(z),function(y) sum(colSums(z[c(x,y),])>0))))
#Jaccard index for the components
coenriched=lapply(1:5,function(x) intersection[[x]]/union[[x]])

#don't forget 1-x to have identical sets together
trees=lapply(coenriched,function(x) hclust(as.dist(1-x)))

#get the groups that are enriched exactly in the same components
groups=lapply(trees,function(x) cutree(x,h=0))
groups=do.call(rbind,lapply(1:5,function(y) 
	data.frame(cbind("subtype"=unique(enriched$BP$subtype)[y],
					 "Description"=names(groups[[y]]),
					 "group"=groups[[y]]))))
write_tsv(groups,"Groups_per_component.tsv")

#copied from http://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Phylogenetic-trees.nb.html
temp=ape::as.phylo(trees$Her2)
ig=ape::as.igraph.phylo(temp,FALSE)
ig=set_edge_attr(ig,"distance",value=temp$edge.length)
createNetworkFromIgraph(ig)
createColumnFilter('junctions', 'id', "^Node\\\\d+$", "REGEX")
junctions<-getSelectedNodes()
setNodeWidthBypass(junctions,1)
setNodeHeightBypass(junctions,1)
setNodeLabelBypass(junctions, "")

#####################CHECK SHARED FUNCTIONS###################3
#ADD GSEA INFO
#files=list.files()
#files=files[grep("gsea",files)]
#gsea=lapply(files,read_tsv)
#names(gsea)=gsub(".gsea","",files)
##match gsea subtype to SGCCA result
#gsea=lapply(gsea,function(x) 
#	x%>%dplyr::select(subtype,Description,NES,p.adjust))
#gsea$KEGG$subtype=gsub("_Normal","",gsea$KEGG$subtype)
#
#j=lapply(gsea,function(x)
#	unique(x$Description[x$p.adjust<0.05]))
##names(j)=names(enriched) check j is named
##  BP KEGG 
## 130   11 
#
##get shared function
#heatmatrix=function(enrichment){
#	#only functions enriched in more than one dataset
#	i=enrichment%>%distinct(subtype,Description)%>%count(Description)%>%
# 	filter(n>1)%>%dplyr::select(Description)%>%unlist
#	enrichment=enrichment[enrichment$Description%in%i,]
#	#count genes & component per function and subtype
#	edges=enrichment%>%dplyr::select(subtype,Description,geneID)%>%
#		group_by(subtype,Description)
#	edges=edges%>%group_map(~length(unique(unlist(strsplit(.x$geneID,
#		"/")))))%>%unlist%>%cbind(group_keys(edges))
#	colnames(edges)[1]="genes"
#	edges1=enrichment%>%group_by(subtype,Description)%>%tally
#	colnames(edges1)[3]="components"
#	edges=merge(edges,edges1,by=c("subtype","Description"))
#return(edges)}
#shared=lapply(enriched,heatmatrix)
##sapply(shared,function(x) length(unique(x$Description)))
##  BP KEGG 
## 233   16 
##BP is too large to plot
#
##for the plot
##paste with heatmatrix output
#temp=lapply(names(shared),function(x) 
#			merge(groups[[x]],shared[[x]],
#				by=c("subtype","Description"),
#				all.y=T))
##ungrouped functions should be the same figure
#temp=lapply(temp,function(x) x%>%group_by(subtype)%>%
#									add_count(group))
#names(temp)=names(shared)
#temp$BP$group[temp$BP$n==1]=0
#temp$KEGG$group[temp$KEGG$n==1]=0
#
#png("KEGGshared.png",width=600)
# ggplot(temp$KEGG)+geom_point(aes(x=subtype,y=Description,
# 	size=components,col=genes,shape=group))+
# scale_color_gradient(low="blue",high="red")+
# theme_light(base_size=16)+scale_size(range=c(3,10))+
# theme(axis.ticks=element_blank())+annotate("text",
# 	y=j$KEGG[j$KEGG%in%shared$KEGG$Description],x="Basal",label="*",
# 	size=7,vjust=.8,hjust=4.3)+coord_cartesian(clip="off")+
# scale_shape_manual(values=c(16:18,15),guide="none")+xlab("")+ylab("")
#dev.off()
#
#BP table is toooo long to plot