#BP matrix for MI
#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)

BP=commandArgs(trailingOnly=TRUE)

transcripts=fread("Her2.transcripts")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),
	filters = 'go',values = BP, mart = mart)
geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP])
transcripts=transcripts[rownames(transcripts)%in%geneSet,]
transcripts$probes=paste("BP",rownames(transcripts),sep='-')
write.table(transcripts,paste("data",BP,sep='/'),sep='\t',
quote=F,row.names=F)
writeLines(geneSet,"bp.txt")