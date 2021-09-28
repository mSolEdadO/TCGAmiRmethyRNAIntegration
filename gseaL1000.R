###################GET L1000 METADATA
library(cTRAP)

meta=loadCMapData("cmapMetadata.txt", type="metadata")
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz
meta <- filterCMapMetadata(meta,perturbationType="trt_cp")
meta <- filterCMapMetadata(meta,cellLine="MDAMB231")
#cell line ~TNB [PMC5665029]

perturbations <- prepareCMapPerturbations(
    metadata=meta,
    zscores="/labs/csbig/multiomics/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
    geneInfo="cmapGeneInfo.txt",compoundInfo="cmapCompoundInfo.txt")

###################GET SETS
library(biomaRt)
#data to translate ensembl ids to symbols
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
	host="http://nov2020.archive.ensembl.org")
#this version has all the genes, the newest dont
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
 			  mart=mart)
#get community data
files=list.files()
files=files[grep("comm",files)]
#merge all in a data.frame
communities=lapply(files,read.csv,header=T)
communities=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(communities[[x]],paste("clust",x,sep='')))))
#translate to hgnc_symbols
communities=communities[order(communities$genes),]
i=table(communities$genes)
communities$symbol=unlist(sapply(1:length(i),function(x) 
	rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(i)[x]],
		i[x])))
communities$name=paste(communities[,3],"community",communities$community,sep='')
#separate by community
sets=lapply(unique(communities$name),function(x) 
	communities$symbol[communities$name==x])
#list elements must be named for fgsea to work
names(sets)=unique(communities$name)

###################ACTUAL ENRICHMENT
library(fgsea)

gseaPerPerturbation<-function(i){
	#make rank
	zscores <- loadCMapZscores(perturbations[ ,i])
	rank=as.numeric(zscores)
	names(rank)=rownames(zscores)
	rank=rank[order(rank,decreasing=T)]
	fgseaRes=fgsea(pathways=sets,stats=rank,minSize=50,nperm=1e5)#nperm recom=10/p.val threshold
#Warning messages:
#There are ties in the preranked stats (0.01% of the list).
#The order of those tied genes will be arbitrary, which may produce unexpected results.
	fgseaRes=cbind(colnames(zscores),fgseaRes[fgseaRes$padj<0.05,])
return(fgseaRes)}

results=do.call(rbind,lapply(1:ncol(perturbations),function(x)
 gseaPerPerturbation(x)))