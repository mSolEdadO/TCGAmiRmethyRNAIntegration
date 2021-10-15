#!/usr/bin/env Rscript
library(cTRAP)

###################GET L1000 METADATA
cells=commandArgs(trailingOnly=TRUE)#like BT20, MDAMB231

#exchange comments for GSE70138
meta=loadCMapData("cmapMetadata.txt", type="metadata")
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz
#meta=loadCMapData("GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt",
#for some reason sig_info files sig_id does NOT match matrix IDs
#        type="metadata")
#meta$cell_id=sapply(strsplit(meta$sig_id,"_"),function(x) x[2])
meta <- filterCMapMetadata(meta,perturbationType="trt_cp")
#trt_cp = Peptides and other biological agents (e.g. cytokine)
meta <- filterCMapMetadata(meta,cellLine=cells)

###################GET SUBSET OF INTEREST
perturbations <- prepareCMapPerturbations(
    metadata=meta,
#change paths to load another GSE
    zscores="/labs/csbig/multiomics/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
#Level 5 (SIG) consists of the replicates, usually 3 per treatment, aggregated into a single differential expression vector derived from the weighted averages of the individual replicates.
#zscores="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
    geneInfo="cmapGeneInfo.txt",compoundInfo="cmapCompoundInfo.txt")
#geneInfo="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",compoundInfo="cmapCompoundInfo.txt")
zscores=loadCMapZscores(perturbations)
write.table(zscores,paste(cells,"mtrx",sep='.'),sep='\t',quote=F)
