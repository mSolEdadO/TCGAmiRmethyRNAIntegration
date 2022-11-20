#!/usr/bin/env Rscript

library(FactoMineR)
library(factoextra)
library(tidyverse)

methy=data.table::fread("methyM.tsv")
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
subtype=read_tsv("subtype.tsv")
subtype=subtype[order(match(subtype$samples,colnames(methy))),]

temp=PCA(t(methy),graph=F)
pdf("methy-pca.pdf")
fviz_pca_ind(temp,col.ind=subtype$subtype,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$treatments_pharmaceutical_treatment_or_therapy,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$treatments_radiation_treatment_or_therapy,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$race,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$ajcc_pathologic_m,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$ajcc_pathologic_n,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$morphology,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$ajcc_pathologic_t,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$prior_malignancy,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$primary_diagnosis,geom.ind="point")
fviz_pca_ind(temp,col.ind=subtype$ajcc_pathologic_stage,geom.ind="point")
dev.off()
