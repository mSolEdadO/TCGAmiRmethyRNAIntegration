library(tftargets)
library(biomaRt)

#############data.frame lists
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
temp=pbapply::pblapply(unique(i),function(x) temp[i==x,])
Neph2012=as.data.frame(do.call(rbind,lapply(temp,function(x)
 cbind(paste(x$V1,collapse=','),x[1,2:3]))))
colnames(Neph2012)=c("tissue","TF","target")

#from 6,175 pubmed articles, with small-scale experiments
#HGNC symbols as targets 
TRRUST=as.data.frame(do.call(rbind,lapply(1:length(TRRUST),
	function(x) cbind(names(TRRUST)[x],TRRUST[[x]]))))

#########################JOIN DBs
ENCODE=cbind("",ENCODE)
colnames(ENCODE)=colnames(Neph2012)
targets=rbind(ENCODE,Neph2012)

#map to same type of ID
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","entrezgene_id",
	"hgnc_symbol"),mart=mart)
#can't keep undefined targets
targets=targets[targets$target%in%myannot$entrezgene_id,]

#collapse dbs by TF-target