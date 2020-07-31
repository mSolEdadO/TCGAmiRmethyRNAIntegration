#!/usr/bin/env Rscript 

library(data.table)
top=fread("Her2.top.tsv")
top=rbind(top,fread("temp"))
top=top[order(top$V2,decreasing=T),]
write.table(top[1:10000,],"Her2.top.tsv",sep='\t',quote=F,row.names=F)

#for x in $(ls|grep 'Her2'|grep 'adj');
#do ../../bin/adj2sif $x>temp;
#Rscript updateTop.R;
#echo $x;
#done