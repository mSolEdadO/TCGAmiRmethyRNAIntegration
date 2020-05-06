library(fgsea)
library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(HTSanalyzeR)

#get databases
GS_KEGG<-KeggGeneSets(species = "Hs")
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    

#get entrez ids of DEGs
mart=useEnsembl("ensembl",
				dataset="hsapiens_gene_ensembl",
				host="http://jan2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene"), 
 			  mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
#focus on transcripts
subDEGs=lapply(de,function(x) x[grep("ENSG",rownames(x)),])
#de is a list with the dfs output of limma's topTreat for different contrasts
temp=lapply(subDEGs,function(x) myannot[myannot$ensembl_gene_id%in%rownames(x),])
temp=lapply(1:4,function(x)
				temp[[x]]$entrezgene[order(match(temp[[x]]$ensembl_gene_id,
					rownames(subDEGs[[x]]),nomatch=0))])

#create needed ranks
ranks=lapply(1:4,function(x) {
				r=subDEGs[[x]]$logFC;
				names(r)=temp[[x]];
				return(r)})
ranks=lapply(ranks,function(x) x[order(x)])
ranks=lapply(ranks,function(x) x[!is.na(names(x))])
#ensembl transcript  mapping to entrez is not one to one
ranks=lapply(ranks,function(x) x[!duplicated(names(x))])
names(ranks)=names(de)

#GSEA
fgseaRes =lapply(ranks,function(x) 
					fgsea(pathways=GS_GO_BP,
						  stats=x, 
		                  minSize=15, 
                  		  nperm=1e5))
sapply(fgseaRes,function(x) summary(x$padj))
#no significant results
fgseaRes =lapply(ranks,function(x) 
					fgsea(pathways=GS_KEGG,
						  stats=x, 
		                  minSize=15, 
                  		  nperm=1e5))
sapply(fgseaRes,function(x) summary(x$padj))
#one significant result: Axon guidance
fgseaRes$luma_normal[fgseaRes$luma_normal$padj<0.05,]
#    pathway         pval       padj         ES       NES nMoreExtreme size
#1: hsa04360 0.0002223155 0.02445471 -0.9266251 -1.851086            8   28
#   leadingEdge
#1:        2049
