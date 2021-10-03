library(cTRAP)

#LOAD metadata
meta <- loadCMapData("cmapMetadata.txt", type="metadata")
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz
meta <- filterCMapMetadata(meta,perturbationType="trt_cp")
drugs <- loadCMapData("cmapCompoundInfo.txt", type="compoundInfo")

#LOAD ENRICHMENT
files=list.files()
files=files[grep("gsea",files)]
data=lapply(files,read.table,header=T)
names(data)=gsub(".gsea","",files)

#ADD METADATA TO ENRICHED SIGNATURES
data=lapply(data,function(x) merge(x,meta[,c(1,3,6:9)],by="sig_id"))
data=lapply(data,function(x) merge(x,drugs[,c(1:6,10)],by="pert_iname"))
data=lapply(data,function(x) x[order(x$NES,decreasing=T),])
#negative NES would be drugs that should not be used, right?
#KEEP ONLY LAUNCHED DRUGS
launched=lapply(data,function(x) x[x$clinical_phase=="Launched",])
write.table(do.call(rbind,launched),"enricheDrugs.tsv",sep='\t',quote=F,row.names=F)
