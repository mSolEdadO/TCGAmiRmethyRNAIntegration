library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(enrichplot)
load("DA.RData")

#translate ensembl ids to ids kegg can use
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene",
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
ranks=lapply(DE.genes,function(x) get_ranks(x,"entrezgene"))
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

###########################
#COMPARE WITH SGCCA OVER-REPRESENTED PATHWAYS
############################
library(tidyverse)
library(ggplot2)

keggenrich=read_tsv("KEGG.enrichment")
png("sgcca_gsea.png")
gseaK=gseaK[gseaK%ID%in%keggenrich$ID,]
gseaK%>%ggplot(aes(subtype,Description,fill=NES,
	alpha=-log(p.adjust)))+
geom_tile()+scale_fill_gradient(low="blue", high="red")+
scale_x_discrete(labels=c("basal_normal"="Basal",
	"her2_normal"="Her2","luma_normal"="LumA","lumb_normal"="LumB"))+
xlab("")+ylab("")+theme(panel.background=element_blank(),
	axis.ticks=element_blank())
dev.off()

gseaK$subtype=gsub("Lumb","LumB",gsub("Luma","LumA",str_to_title(gsub("_normal","",gseaK$subtype))))
gseaK%>%filter(p.adjust<0.01)%>%select(subtype,ID,NES)%>%merge(keggenrich,by=c("subtype","ID"))
