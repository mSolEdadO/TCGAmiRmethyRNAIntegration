library(data.table)
library(ggplot2)
#load sif files obtained with adj2sif
files=list.files("parallel-aracne",full.names=T)
files=files[grep("sif",files)]
sif=lapply(files,fread)
names(sif)=gsub("parallel-aracne/","",gsub(".miR.sif","",files))
#########################TOP MI interactions#########################
sif=lapply(sif,function(x) x[order(x$V2,decreasing=T),])
i=lapply(sif,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
j=unique(i$Basal)
top=lapply(1:5,function(x) 
	do.call(rbind,lapply(j,function(y) sif[which(i[[x]]==y)[1:10000],])))
i=lapply(top,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
