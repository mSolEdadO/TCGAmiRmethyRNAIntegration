#!/usr/bin/env Rscript 
library(data.table)

files=commandArgs(trailingOnly=TRUE)

ori=fread(files[1])
data=fread(files[2])
#print(nrow(data))
#data$index=1:nrow(data)
#get indexes
o=paste(ori$V1,ori$V3)
o1=paste(ori$V3,ori$V1)
i=paste(data$V1,data$V2)

#original interactions in this iteration
found=rbind(data[i%in%o,],data[i%in%o1,])
i=paste(found$V1,found$V2)
temp=found$V2[!i%in%o]
found$V2[!i%in%o]=found$V1[!i%in%o]
found$V1[!i%in%o]=temp
i=paste(found$V1,found$V2)
found=found[order(match(i,o)),]
final=gsub(".sif",".sort",files[2])
write.table(found,final,sep='\t',quote=F,row.names=F)