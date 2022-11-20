#!/usr/bin/env Rscript
######BEFORE concatena.R coz
#var distributions should be similar if the planned 
#approach relies on correlation networks [Tarazona2020]

library(data.table)#1.14.2
library(FactoMineR)#2.4
library(factoextra)#1.0.7
library(ggplot2)#3.3.5

omic = commandArgs(trailingOnly=TRUE)

data=fread(paste(omic,"tsv",sep='.'))
data=t(as.matrix(data[,2:ncol(data)],rownames=data$V1))
#factominer eats individuals in rows & variables in columns

####THE NORMALIZATION
pcs=PCA(data, scale.unit = TRUE, ncp = 2, graph = F)
#print("eigenvalue > 1 indicates the PC accounts for more variance than original variables")
#print(sum(pcs$eig[,1]>1))
w=pcs$svd$vs[1]
normi=t(data)/w
write.table(normi,paste(omic,"eigenNormi",sep='.'),sep='\t',quote=F)


