library(fgsea)
library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(HTSanalyzeR)
load("DA.RData")#DA.R output
#get databases
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    

#get entrez ids of DEGs
mart=useEnsembl("ensembl",
				dataset="hsapiens_gene_ensembl",
				host="http://jan2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene"), 
 			  mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
temp=lapply(DE.genes,function(x) myannot[myannot$ensembl_gene_id%in%rownames(x),])
temp=lapply(1:4,function(x)
				temp[[x]]$entrezgene[order(match(temp[[x]]$ensembl_gene_id,
					rownames(DE.genes[[x]]),nomatch=0))])

#create needed ranks
ranks=lapply(1:4,function(x) {
				r=DE.genes[[x]]$t;
#as suggested in https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
				names(r)=temp[[x]];
				return(r)})
ranks=lapply(ranks,function(x) x[order(x)])
ranks=lapply(ranks,function(x) x[!is.na(names(x))])
#ensembl transcript  mapping to entrez is not one to one
ranks=lapply(ranks,function(x) x[!duplicated(names(x))])
names(ranks)=names(DE.genes)

#GSEA
fgseaRes =lapply(ranks,function(x) 
					fgsea(pathways=GS_GO_BP,
						  stats=x, 
		                  minSize=15, 
                  		  nperm=1e5))#nperm recom=10/p.val threshold

fgseaRes=lapply(fgseaRes,function(x) x[x$padj<0.01,])
sapply(fgseaRes,nrow)
#basal_normal  her2_normal  luma_normal  lumb_normal       normal 
#         182          165           95          217          270 
