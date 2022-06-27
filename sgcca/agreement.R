library(tidyverse)
library(irr)

files=list.files()
subsa=lapply(files,read_tsv)
final=read_tsv("../LumB.selected")

subsa=lapply(subsa,function(x) x%>%count(feature))
subsa=reduce(subsa,full_join,by="feature")
colnames(subsa)[2:101]=paste("r",1:100,sep="")
final=final%>%count(variable)
colnames(final)[1]="feature"
subsa=merge(final,subsa,by="feature",all.x=T)

temp=apply(subsa[,2:101],1,table)
write_tsv(subsa,"LumB.subsample.selected")
