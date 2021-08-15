#BP matrix for MI
#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)

BP=commandArgs(trailingOnly=TRUE)

transcripts=fread("Her2.transcripts")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
if(length(mart)==0){stop("biomaRt didn't work")}
geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),
	filters = 'go',values = BP, mart = mart)
geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP])
transcripts=transcripts[transcripts$probes%in%geneSet,]
transcripts$probes=paste("BP",transcripts$probes,sep='-')
write.table(transcripts,file=BP,sep='\t',quote=F,row.names=F)
writeLines(geneSet,"../bp.txt")