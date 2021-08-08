library(tftargets)
library(biomaRt)

#############################TABLE FROM LISTS
#entrez ID for targets
#based on ChIP-seq 
ENCODE=as.data.frame(do.call(rbind,lapply(1:length(ENCODE),
	function(x) cbind(names(ENCODE)[x],ENCODE[[x]]))))

#discovered by DNaseI footprinting and TF seq 
temp=lapply(Neph2012,function(x) x[lapply(x,length)>0])
#per tissue
temp=lapply(1:length(temp),function(x) do.call(rbind,
	lapply(1:length(temp[[x]]),function(y) cbind(names(temp[[x]])[y],
		temp[[x]][[y]]))))
names(temp)=names(Neph2012)
temp=as.data.frame(do.call(rbind,lapply(1:length(temp),function(x) 
	cbind(names(temp)[x],temp[[x]]))))
#collapse tissues by TF-target
i=paste(temp$V2,temp$V3)
temp=lapply(unique(i),function(x) temp[i==x,])
Neph2012=as.data.frame(do.call(rbind,lapply(temp,function(x)
 cbind(paste(x$V1,collapse=','),x[1,2:3]))))
colnames(Neph2012)=c("tissue","TF","target")

#small-scale experiments from pubmed
#HGNC symbols for targets 
TRRUST=as.data.frame(do.call(rbind,lapply(1:length(TRRUST),
	function(x) cbind(names(TRRUST)[x],TRRUST[[x]]))))

#########################JOIN DBs
ENCODE=cbind("",ENCODE,"ENCODE")
Neph2012$DB="Neph2012"
colnames(ENCODE)=colnames(Neph2012)
targets=rbind(ENCODE,Neph2012)

#map to same type of ID
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"hgnc_symbol"),mart=mart)
myannot=myannot[!is.na(myannot$entrezgene_id),]

#can't keep undefined targets
sum(!targets$target%in%myannot$entrezgene_id)
#[1] 580
targets=targets[targets$target%in%myannot$entrezgene_id,]
targets=targets[order(targets$target),]
i=table(targets$target)
#keep all ensembl ids that match to the same entrez id
targets$target=unlist(sapply(1:length(i),function(x)
	rep(paste(unique(myannot$ensembl_gene_id[myannot$entrezgene_id%in%
		names(i)[x]]),collapse=','),i[x])))
sum(!TRRUST$V2%in%myannot$hgnc_symbol)
#[1] 82
TRRUST=TRRUST[TRRUST$V2%in%myannot$hgnc_symbol,]
TRRUST=TRRUST[order(TRRUST$V2),]
i=table(TRRUST$V2)
#keep all ensembl ids that match to the same entrez id
TRRUST$V2=unlist(sapply(1:length(i),function(x)
	rep(paste(unique(myannot$ensembl_gene_id[myannot$hgnc_symbol%in%
		names(i)[x]]),collapse=','),i[x])))

#1 line per TF-target
TRRUST=cbind("",TRRUST,"TRRUST")
colnames(TRRUST)=colnames(targets)
targets=rbind(unique(TRRUST),targets)#TRRUS has repated rows
#the same ensembl_id matches several entrez_id??? 
i=paste(targets$TF,targets$target)
temp=unique(targets[i%in%i[duplicated(i)],])
targets=targets[!i%in%i[duplicated(i)],]
i=paste(temp$TF,temp$target)
temp=lapply(unique(i),function(x) temp[i==x,])
temp1=as.data.frame(do.call(rbind,lapply(temp,function(x)
  cbind(x[nrow(x),1:3],paste(unique(x$DB),collapse=',')))))#only Neph2012 has tissue info & is always the last row in the list
colnames(temp)=colnames(targets)
targets=rbind(targets,temp)

#can't use TFs with no annotation
sum(!targets$TF%in%myannot$hgnc_symbol)
#[1] 918 #may change with ensembl version
targets=targets[targets$TF%in%myannot$hgnc_symbol,]
targets=targets[order(targets$TF),]
i=table(targets$TF)
targets$TF=unlist(sapply(1:length(i),function(x) 
	rep(paste(unique(myannot$ensembl_gene_id[myannot$hgnc_symbol%in%
		names(i)[x]]),collapse=","),i[x])))
write.table(targets,"TFtargets.tsv",sep='\t',quote=F,row.names=F)