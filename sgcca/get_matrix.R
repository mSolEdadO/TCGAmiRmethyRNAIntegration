#single dataframe with the 4 subtypes & the normal tissue
sets=do.call(rbind,lapply(1:5,function(x) 
	cbind(subtype=names(sets)[x],
		sets[[x]][,c("component","variable")])))

#get functions enriched in the 4 subtypes & the normal tissue
shared=sapply(functions,function(x) 
	names(which(table(unlist(x))==5)))
enriched=lapply(1:2,function(x) 
	enriched[[x]][enriched[[x]]$ID%in%shared[[x]],])
names(enriched)=names(shared)
#drop child terms 
todrop=c("GO:0007156",#GO:0007156 is a GO:0098742
	"GO:0072009",#is a GO:0072006 & GO:0072073                              
	"GO:0072080")#is a GO:0061326 & GO:0072009
enriched$BP=enriched$BP%>%filter(!ID%in%todrop)

#get features for selected functions
features=lapply(enriched,function(x) 
	x%>%dplyr::select(subtype,component,ID)%>%
	merge(sets,by=c("subtype","component")))
