#!/usr/bin/env Rscript
library(tidyverse)

########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
fun=args[1]
subty=args[2]

#get the components linked to the function
if(length(grep("GO",fun))>0){
 enrich=read_tsv("BP.enrichment")
}else{
 enrich=read_tsv("KEGG.enrichment")	
}
comp=enrich%>%filter(ID==fun&subtype==subty)%>%
			  dplyr::select(component)%>%unlist

#get the features selected in those components
selected=read_tsv(paste(subty,"stable",sep='.'))
features=selected%>%filter(component%in%comp)%>%
		 distinct(variable)%>%unlist
#print(paste("Function has",length(features),"features associated",sep=' '))

#get the data
data=data.table::fread(paste(subty,"eigeNormi",sep='.'))
data=data[data$V1%in%features,]
write_tsv(data,paste(fun,subty,"mtrx",sep='.'))

