#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
BP=args[1]
type=args[2]

targets=readLine("bp.txt")
sub=readLine("condor.seed")

sub=gsub("GO",BP,sub)
sub1=sapply(bp,function(x) gsub("ENSG",x))
sub2=gsub(paste(BP,"txt",sep='.'),
	paste(BP,"cpgs.txt",sep='-'),
	sub1)
sub2=gsub("-p 0.001","-p 0.000001",sub1)
writeLines(c(sub1,sub2),paste(BP,"sub",sep='.'))