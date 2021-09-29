#!/usr/bin/env Rscript
library(cTRAP)

###################GET L1000 METADATA
cells=commandArgs(trailingOnly=TRUE)

meta=loadCMapData("cmapMetadata.txt", type="metadata")
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz
meta <- filterCMapMetadata(meta,perturbationType="trt_cp")
meta <- filterCMapMetadata(meta,cellLine=cells)

###################GET SUBSET OF INTEREST
perturbations <- prepareCMapPerturbations(
    metadata=meta,
    zscores="/labs/csbig/multiomics/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
#Level 5 (SIG) consists of the replicates, usually 3 per treatment, aggregated into a single differential expression vector derived from the weighted averages of the individual replicates.
    geneInfo="cmapGeneInfo.txt",compoundInfo="cmapCompoundInfo.txt")
#change zscores path to load another GSE
zscores=loadCMapZscores(perturbations)
write.table(zscores,paste(cells,"mtrx",sep='.'),sep='\t',quote=F)
