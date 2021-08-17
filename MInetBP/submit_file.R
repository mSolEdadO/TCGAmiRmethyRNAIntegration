#!/usr/bin/env Rscript

BP=commandArgs(trailingOnly=TRUE)

targets=readLines("bp.txt")
sub=readLines("condor.seed")

sub=gsub("GO",BP,sub)
targets=paste("BP",targets,sep='-')
sub1=sapply(targets,function(x) gsub("ENSG",x,sub))
sub2=gsub(BP,paste(BP,"cpgs",sep='-'),sub1)
sub2=gsub("-p 0.001","-p 0.000001",sub2)
writeLines(c(sub1,sub2),paste(BP,"sub",sep='.'))