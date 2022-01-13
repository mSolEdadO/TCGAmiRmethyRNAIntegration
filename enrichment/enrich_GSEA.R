library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(enrichplot)
library(tidyverse)
load("DA.RData")

#translate ensembl ids to ids kegg can use
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"entrezgene_accession"), mart=mart)
DE.genes=lapply(DE.genes,function(x) 
	cbind(ensembl_gene_id=rownames(x),x))
DE.genes=lapply(DE.genes,function(x) 
	merge(x,myannot,by="ensembl_gene_id"))
#ranked genelist
get_ranks=function(x,ID){ 
	r=as.numeric(x$t);
	names(r)=unlist(x[ID]);
	r=r[order(r,decreasing=T)];
	r=r[!is.na(names(r))];
	r=r[!duplicated(names(r))];
	return(r)}
ranks=lapply(DE.genes,function(x) get_ranks(x,"entrezgene_id"))
#if ensembl_gene_id R crashes

#GO GSEA
gseaGO <-lapply(ranks,function(x) 
	gseGO(geneList     = x,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
           # keyType       = 'ENSEMBL',
            pvalueCutoff = 1,#to recover all sgcca over-represented paths 
			pAdjustMethod="fdr"))
gseaGO=do.call(rbind,lapply(1:4,function(x) 
	cbind(subtype=names(gseaK)[x],as.data.frame(gseaK[[x]]))))

#KEGG GSEA
gseaK=lapply(ranks,function(x) 
	gseKEGG(geneList=x,
		organism='hsa',
		pvalueCutoff = 1,#to recover all sgcca over-represented paths 
		pAdjustMethod="fdr"))
gseaK=do.call(rbind,lapply(1:4,function(x) 
	cbind(subtype=names(gseaK)[x],as.data.frame(gseaK[[x]]))))

write_tsv(gseaK,"KEGG.gsea")