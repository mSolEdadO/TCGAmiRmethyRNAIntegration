library(clusterProfiler)
library(org.Hs.eg.db)
#library(biomaRt)
library(enrichplot)
library(tidyverse)


#translate ensembl ids to ids kegg can use
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=105)
#myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
#	"entrezgene_accession"), mart=mart)
DE.genes=read_tsv("DE.genes.tsv")
#DE.genes=merge(DE.genes,myannot,by="ensembl_gene_id")
#choosen over bitr, which fails to map 0.14% of input gene IDs

#ranked genelist
get_ranks=function(x,ID){ 
	r=as.numeric(x$logFC);
	names(r)=unlist(x[ID]);
	r=r[order(r,decreasing=T)];
	r=r[!is.na(names(r))];
	r=r[!duplicated(names(r))];
	return(r)}
ranks=lapply(unique(DE.genes$contrast),function(x) 
	get_ranks(DE.genes[DE.genes$contrast==x,],"entrezgene_id"))
#if ensembl_gene_id R crashes
names(ranks)=unique(DE.genes$contrast)

#GO GSEA
#gseaGO=list()
#gseaGO$LumB=gseGO(geneList     = ranks$LumB_Normal,
#            OrgDb        = org.Hs.eg.db,
#            ont          = "BP",
#           # keyType       = 'ENSEMBL',
#            pvalueCutoff = 1,#to recover all sgcca over-represented paths 
#			pAdjustMethod="fdr")
#
gseaGO <-lapply(ranks,function(x) 
	gseGO(geneList     = x,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
           # keyType       = 'ENSEMBL',
            pvalueCutoff = 1,#to recover all sgcca over-represented paths 
			pAdjustMethod="fdr"))
gseaGO=do.call(rbind,lapply(1:4,function(x) 
	cbind(subtype=names(gseaGO)[x],as.data.frame(gseaGO[[x]]))))
write_tsv(gseaGO,"BP.gsea")

#KEGG GSEA
gseaK=lapply(ranks,function(x) 
	gseKEGG(geneList=x,
		organism='hsa',
		pvalueCutoff = 1,#to recover all sgcca over-represented paths 
		pAdjustMethod="fdr"))
gseaK=do.call(rbind,lapply(1:4,function(x) 
	cbind(subtype=names(gseaK)[x],as.data.frame(gseaK[[x]]))))

write_tsv(gseaK,"KEGG.gsea")