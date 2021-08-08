library(fgsea)
library(biomaRt)
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(HTSanalyzeR)
#load("DA.RData")#DA.R output: DE.genes.tsv
#get databases
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))                    

#get entrez ids of DEGs
mart=useEnsembl("ensembl",
				dataset="hsapiens_gene_ensembl",
				host="http://jan2019.archive.ensembl.org")
#drop host for latest realease
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene"),
#entrezgene → entrezgene_id in latest ensembl 
 			  mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
myannot=myannot[myannot$ensembl_gene_id%in%DE$gene,]
temp=lapply(unique(DE$contrast),function(x) DE[DE$contrast==x,])
temp=lapply(temp,function(x) x[order(x$gene),])
#drop DE genes with no entrez id
i=unique(DE$gene[!DE$gene%in%myannot$ensembl_gene_id])
temp=lapply(temp,function(x) x[!x$gene%in%i,])
myannot=myannot[order(match(myannot$ensembl_gene_id,temp[[1]]$gene)),]

#create needed ranks
ranks=lapply(1:4,function(x) {
				r=temp[[x]]$t;
#as suggested in https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
				names(r)=myannot$entrezgene_id;
				return(r)})
ranks=lapply(ranks,function(x) x[order(x)])
ranks=lapply(ranks,function(x) x[!is.na(names(x))])
#ensembl transcript  mapping to entrez is not one to one
ranks=lapply(ranks,function(x) x[!duplicated(names(x))])
names(ranks)=unique(DE$contrast)

#GSEA
fgseaRes =lapply(ranks,function(x) 
					fgsea(pathways=GS_GO_BP,
						  stats=x, 
		                  minSize=15, 
                  		  nperm=1e5))#nperm recom=10/p.val threshold

fgseaRes=lapply(fgseaRes,function(x) x[x$padj<0.05,])
sapply(fgseaRes,nrow)
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#          99           76           53          108 
#keep only parent ones
leading=lapply(1:4,function(x) 
	collapsePathways(fgseaRes[[x]],GS_GO_BP,ranks[[x]]))
fgseaRes=lapply(1:4,function(x) 
	fgseaRes[[x]][fgseaRes[[x]]$pathway%in%leading[[x]]$mainPathways,])
sapply(fgseaRes,nrow)
#[1] 67 47 32 63
fgseaRes=as.data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(unique(DE$contrast)[x],fgseaRes[[x]]))))
fgseaRes$pathName=Term(fgseaRes$pathway)
colnames(fgseaRes)[1]="contrast"
write.table(apply(fgseaRes,2,as.character),"fgseaRes.tsv",
	sep='\t',quote=F,row.names=F)

#plot
library(igraph)
library(gplots)

g=graph.data.frame(fgseaRes,c(1,9)],directed=F)
E(g)$w=as.numeric(as.character(fgseaRes$NES))
A=get.adjacency(g,attr="w")
A=as.matrix(A[unique(fgseaRes$pathName),unique(fgseaRes$contrast)])
pdf("NES.pdf",height=15)
heatmap.2(A,trace='n',scale='n',col=greenred(10),dendrogram='r',
	Colv=F,key.title=F,key.xlab="NES",keysize=0.25,
	lmat = rbind(c(0,4),c(2,1),c(0,3)),lhei=c(0.5,5,0.2),
	lwid=c(1,6),srtCol=0,adjCol = c(0.5,1),
	labCol=gsub("_Normal","",colnames(A)),margins=c(2,25))
dev.off()