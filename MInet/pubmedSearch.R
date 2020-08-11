library(pbapply)
library(data.table)
library(rentrez)
set_entrez_key("49b3079321d573aaa12522e38a1b31d38e08")#ncbi account for dopreto to submit 10 queries per second

#########################CpGs#########################
methy=fread("MapMethy.tsv")
#keep only CpGs per BP per subtype
cpgs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="c"]))
#build queries
queries=unique(unlist(lapply(cpgs,function(x) lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	methy$UCSC_RefGene_Name[methy$IlmnID%in%x[[y]]],
	    	"(CpG methylation OR DNA methylation)",
	    	sep=" AND ")))))
length(queries)

#search db
comention=pblapply(queries,function(x) {
		  request=entrez_search(db = "pubmed", term = y);
		  Sys.sleep(0.01);
		  return(request)})
#get succesful queries
tfs=unlist(comention[2,unlist(comention[2,])>0])

#########################TFs#########################
#all over again
tfs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="E"]))
queries=unique(unlist(lapply(tfs,function(x) lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	myannot$hgnc_symbol[myannot$ensembl_gene_id%in%x[[y]]],
	    	sep=" AND ")))))
length(queries)
#[1] 15683
comention=pbsapply(queries,function(x) {
		  request=entrez_search(db = "pubmed", term = x);
		  Sys.sleep(0.01);
		  return(request)})
tfs=unlist(comention[2,unlist(comention[2,])>0])
length(tfs)

#########################miRNAs#########################
#all over again
mirs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="h"]))
queries=unique(unlist(lapply(mirs,function(x) lapply(1:length(x),function(y) 
	    paste(names(x)[y],x[[y]],sep=" AND ")))))
length(queries)
#[1] 22935
comention=pbsapply(queries,function(x) {
		  request=entrez_search(db = "pubmed", term = x);
		  Sys.sleep(0.01);
		  return(request)})
mirs=unlist(comention[2,unlist(comention[2,])>0])
length(mirs)
#[1] 2165
