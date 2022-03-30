#!/usr/bin/env Rscript
library(tidyverse)

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
fun=args[1]
subty=args[2]

#get the components linked to the function
if(grep("GO",fun)){
 enrich=read_tsv("BP-allFeatures.enrichment")
}else{
 enrich=read_tsv("KEGG-allFeatures.enrichment")	
}
comp=enrich%>%filter(ID==fun&subtype==subty)%>%
			  dplyr::select(component)%>%unlist

#get the features selected in those components
selected=read_tsv(paste(subty,"selected",sep='.'))
features=selected%>%filter(component==unlist(comp))%>%
		 distinct(variable)%>%unlist
#print(paste("Function has",length(features),"features associated",sep=' '))

#get the data
data=data.table::fread(paste(subty,"eigeNormi",sep='.'))
data=data[data$V1%in%features,]
write_tsv(data,paste(fun,subty,"mtrx",sep='.'),col_names=F)#needed for puma

