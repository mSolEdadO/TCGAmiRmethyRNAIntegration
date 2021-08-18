#!/usr/bin/env Rscript

BP=commandArgs(trailingOnly=TRUE)#GO:XXX of interest

targets=readLines("bp.txt")#transcripts linked to BP
sub=readLines("condor.seed")#executable=/usr/bin/aracne	arguments= -i data/GO.txt -h ENSG -p 0.000001	output=GO.ENSG.output	error=GO.ENSG.error	log=GO.ENSG.log	queue
#-p values attempted to keep all BP nodes, check for larger BP

#update to BP
sub=gsub("GO",BP,sub)
#to identify transcripts linked to BP
targets=paste("BP",targets,sep='-')
#get condor lines for each transcript in the BP
sub1=sapply(targets,function(x) gsub("ENSG",x,sub))
#get condor lines for the matrix with CpGs
sub2=gsub(BP,paste(BP,"cpgs",sep='-'),sub1)
#more restrictive -p so u get less edges
sub2=gsub("-p 0.001","-p 0.000001",sub2)
writeLines(c(sub1,sub2),paste(BP,"sub",sep='.'))