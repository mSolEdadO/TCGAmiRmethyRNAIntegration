library(pbapply)
library(data.table)
library(rentrez)
set_entrez_key("49b3079321d573aaa12522e38a1b31d38e08")#ncbi account for dopreto to submit 10 queries per second

#########################CpGs#########################
methy=fread("MapMethy.tsv")
#keep only CpGs per BP per subtype
cpgs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="c"]))
#build queries
queries=lapply(cpgs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	methy$UCSC_RefGene_Name[methy$IlmnID%in%x[[y]]],
	    	"(CpG methylation OR DNA methylation)",
	    	sep=" AND "))))
#search db
comention=pblapply(queries,function(x) sapply(x,function(y) {
		  request=entrez_search(db = "pubmed", term = y);
		  Sys.sleep(0.1);
		  return(request)}))
#get index of succesful queries
count=lapply(comention,function(x) apply(x,2,function(y) y$count))
#keep succesful queries
cpgs=lapply(1:5,function(x) cbind(queries[[x]][count[[x]]>0],
							count[[x]][count[[x]]>0]))

#########################TFs#########################
#all over again
tfs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="E"]))
queries=lapply(tfs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	myannot$hgnc_symbol[myannot$ensembl_gene_idx[[y]]],
	    	sep=" AND "))))
comention=pblapply(queries,function(x) sapply(x,function(y) {
		  request=entrez_search(db = "pubmed", term = y);
		  Sys.sleep(0.1);
		  return(request)}))
count=lapply(comention,function(x) apply(x,2,function(y) y$count))
tfs=lapply(1:5,function(x) cbind(queries[[x]][count[[x]]>0],
							count[[x]][count[[x]]>0]))

#########################miRNAs#########################
#all over again
mirs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="h"]))
queries=lapply(mirs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],x[[y]],sep=" AND "))))
#attemp to accelerate
allQ=unique(unlist(queries))
comention=pbsapply(allQ,function(x) {
		  request=entrez_search(db = "pubmed", term = x);
		  Sys.sleep(0.1);
		  return(request)})
count=lapply(comention,function(x) apply(x,2,function(y) y$count))
mirs=lapply(1:5,function(x) cbind(queries[[x]][count[[x]]>0],
							count[[x]][count[[x]]>0]))
