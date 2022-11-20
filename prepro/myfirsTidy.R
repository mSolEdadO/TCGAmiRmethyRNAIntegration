library(tidyverse)

files=list.files()
files=files[grep("ecart",files)]
data=lapply(files,read_tsv)
names(data)=gsub(".ecart","",files)
#instead of a list use tidy groups
data=do.call(rbind,lapply(1:4,function(x) 
		cbind(data[[x]],names(data)[x])))
data%>%group_by(subtype)

#unifor methylation is not of interest
cgs=data%>%filter(group==1)
cgs%>%filter(ecart.type==0)
