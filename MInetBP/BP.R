#BP matrix for MI
#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)

BP=commandArgs(trailingOnly=TRUE)#GO:XXX of interest
transcripts=fread("Her2.transcripts")#update per subtype

#get transcripts linked to BP
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
#biomaRt version fails a lot, so it dies when output aint achieved
if(length(mart)==0){stop("biomaRt didn't work")}
geneSet<- getBM(attributes=c('ensembl_gene_id', 'go_id'),
	filters = 'go',values = BP, mart = mart)
geneSet=unique(geneSet$ensembl_gene_id[geneSet$go_id==BP])
#matrix with the transcripts of interest
transcripts=transcripts[transcripts$probes%in%geneSet,]
transcripts$probes=paste("BP",transcripts$probes,sep='-')
write.table(transcripts,file=BP,sep='\t',quote=F,row.names=F)
writeLines(geneSet,"../bp.txt")