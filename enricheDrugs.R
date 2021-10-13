library(cTRAP)

#LOAD metadata
meta <- loadCMapData("cmapMetadata.txt", type="metadata")
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz
#keep  IDs
ids <- meta[,c("sig_id","pert_iname")]
meta <- loadCMapData("GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt",
		type="metadata")
ids=rbind(ids,meta[,c("sig_id","pert_iname")])
drugs <- loadCMapData("cmapCompoundInfo.txt", type="compoundInfo")

#LOAD ENRICHMENT
files<-list.files()
files<-files[grep("gsea",files)]
data<-lapply(files,read.table,header=T)
names(data)<-gsub(".gsea","",files)

#ADD METADATA TO ENRICHED SIGNATURES
data=lapply(data,function(x) merge(x,ids,by="sig_id"))
data=lapply(data,function(x) merge(x,drugs,by="pert_iname"))
data=lapply(data,function(x) x[order(x$NES,decreasing=T),])

#KEEP ONLY LAUNCHED DRUGS
launched=lapply(data,function(x) x[x$clinical_phase=="Launched",])
#separate per series and tissue
write.table(do.call(rbind,launched[grep("GSE",names(launched))]),
	"BreastenricheDrugs_GSE70138.tsv",sep='\t',quote=F,row.names=F)

cells=c("BT20","MDAMB231","HS578T","MCF10A")
launched=launched[grep("GSE",names(launched),invert=T)]
write.table(do.call(rbind,launched[names(launched)%in%cells]),
	"BreastenricheDrugs",sep='\t',quote=F,row.names=F)
#lung results
write.table(do.call(rbind,launched[!names(launched)%in%cells]),
	"LungEnricheDrugs",sep='\t',quote=F,row.names=F)


########################ARE COMMUNITIES THE SAME ACROSS NETS???
#library(biomaRt)
#library(HTSanalyzeR)
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://apr2019.archive.ensembl.org")
#myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","mirbase_id","external_gene_name","entrezgene"), mart=mart)
#used ensembl ids since communities original ids are of this type
#universo=as.character(myannot$ensembl_gene_id)
#communities=read.table("communities.tsv",sep='\t',header=T)
#communities$group=gsub("community[0-9]+","",communities$name,perl=T)
#sets=lapply(unique(communities$group),function(x) communities[communities$group==x,])
#sets=lapply(sets,function(x) 
#	lapply(unique(x$name),function(y) x$genes[x$name==y]))
#sameSet=temp=lapply(1:3,function(x) lapply((x+1):4,function(z) 
#	do.call(rbind,lapply(sets[[x]],function(y)
#	 multiHyperGeoTest(collectionOfGeneSets=sets[[z]],
#	 	universe=universo,hits=as.character(y),minGeneSetSize=15,
#	 	pAdjustMethod="fdr")))))
#modules contain the same genes
#sapply(sameSet,function(x) sapply(x,function(y) sum(y[,7]<0.05)))
#[[1]]
#[1] 23 24 24
#[[2]]
#[1] 25 23
#[[3]]
#[1] 22
