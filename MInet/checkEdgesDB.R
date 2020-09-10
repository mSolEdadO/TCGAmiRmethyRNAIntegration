library(data.table)
library(pbapply)
library(multiMiR)

#filter out function-transcript interactions
i=sapply(regus,function(x) names(x)[1])
i=sapply(1:5,function(x) which(BPenriched[[x]]$V1==i[x])[1]-1)
temp=lapply(1:5,function(x) BPenriched[[x]][1:i[x],])
#separate each type of interactions
i=lapply(temp,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
cpgs=lapply(1:5,function(x) temp[[x]][i[[x]]=="E c",])
tfs=lapply(1:5,function(x) temp[[x]][i[[x]]=="E E",])
#not all transcripts are TFs
tfs=lapply(tfs,function(x) x[x$V1%in%unlist(regus)|x$V3%in%unlist(regus),])
mirs=lapply(1:5,function(x) temp[[x]][i[[x]]=="h E",])

#########################CpGs#########################
#Illumina file
methy=fread("HumanMethylation450_15017482_v1-2.csv",sep=',',header=T,fill=T,skip=7)
#only keep those with position and gene info
methy=methy[!is.na(methy$MAPINFO),c(1,12,13,22,24)]
methy=methy[methy$UCSC_RefGene_Name!="",]
#1 line per gene
methy=data.frame(do.call(rbind,apply(methy,1,function(x) 
	cbind(paste(x[1:3],collapse='+'),
		  unlist(strsplit(x[4],";")),
		  unlist(strsplit(x[5],";"))))))
methy=cbind(do.call(rbind,strsplit(as.character(methy$X1),"+",fixed=T)),methy[,2:3])
colnames(methy)=c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")
write.table(methy,"MapMethy.tsv",sep='\t',quote=F,row.names=F)

#recover methy lines with cpgs in enriched sifs
true.cpgs=pblapply(cpgs,function(y) do.call(rbind,apply(y,1,function(x) 
	methy[methy$UCSC_RefGene_Name==myannot$hgnc_symbol[myannot$ensembl_gene_id==x[1]]&
		  methy$IlmnID==x[3],])))
#count unique cpg-gene relations
sapply(true.cpgs,function(x) nrow(unique(x[,c(1,4)])))
#[1] 554  88 536 708  28
#same chr CpG-genes
temp=unique(methy[,1:2])
true.cpgs=pbsapply(cpgs,function(y) sum(apply(y,1,function(x) 
	myannot$chromosome_name[myannot$ensembl_gene_id==x[1]]==temp$CHR[temp$IlmnID==x[3]])))
sapply(true.cpgs,function(x) sum(unlist(x)))
#[1] 1269  435  938 2145  235

#########################TFs#########################
tftargets=fread("TFtargets.tsv")#from https://github.com/slowkow/tftargets

#TF could be both in V1 or V3
side1=pblapply(tfs,function(x) do.call(rbind,apply(x,1,function(y) 
	tftargets[tftargets$V3==myannot$hgnc_symbol[myannot$ensembl_gene_id==y[1]]&
	tftargets$V4==myannot$hgnc_symbol[myannot$ensembl_gene_id==y[3]],])))
side2=pblapply(tfs,function(x) do.call(rbind,apply(x,1,function(y) 
	tftargets[tftargets$V4==myannot$hgnc_symbol[myannot$ensembl_gene_id==y[1]]&
	tftargets$V3==myannot$hgnc_symbol[myannot$ensembl_gene_id==y[3]],])))
true.tfs=lapply(1:5,function(x) rbind(side1[[x]],side2[[x]]))
#TRRUST interactions are supported by focused articles manually curated
sapply(true.tfs,function(x) nrow(unique(x[x$V1=="TRRUST",3:4])))
#[1]  5  2  5  1 14
#the rest are not validated but support in some way the interaction
sapply(true.tfs,function(x) nrow(unique(x[x$V1!="TRRUST",3:4])))
#[1] 127 121 135 133 282

#########################miRNAs#########################
#I have precursor IDs but multimir uses mature IDs
mirIDs=fread("miR.ids.map.tsv",skip=1)#from #ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.xls.gz
#some of my mirnas are not in the table coz they're not in release 22
sum(!unique(unlist(sapply(mirs,function(x) x$V1)))%in%mirIDs$precursor)
#[1] "hsa-mir-4461" "hsa-mir-3607" "hsa-mir-3653"

#get all interactions of the all miRNAs
miRtargets=get_multimir(mirna=mirIDs$mature,summary=F,table="all",legacy.out=F)
miRtargetsV=select(miRtargets,keys="validated",columns=columns(miRtargets),keytype="type")
miRtargetsP=select(miRtargets,keys="predicted",columns=columns(miRtargets),keytype="type")
miRtargets=rbind(miRtargetsV,miRtargetsP)

true.mirs=pblapply(mirs,function(x) do.call(rbind,apply(x,1,function(y)
 miRtargets[miRtargets$mature_mirna_id==mirIDs$mature[mirIDs$precursor==y[1]]&
 miRtargets$target_ensembl==y[3],])))
temp=lapply(true.mirs,function(x) unique(x[,c(3,6,10)]))
sapply(temp,function(x) table(x[,3]))
#          [,1] [,2] [,3] [,4] [,5]
#predicted  236  304  177  339  445
#validated  167  226  111  201  284
temp1=lapply(temp,function(x) x[duplicated(x[,1:2]),])
sapply(temp,function(x) table(x[,3]))[2,]-sapply(1:5,function(x) 
	sum(apply(temp1[[x]],1,function(y) 
		temp[[x]][temp[[x]][,1]==y[1]&temp[[x]][,2]==y[2],3])=="validated"))
#[1] 115 159  68 143 198
