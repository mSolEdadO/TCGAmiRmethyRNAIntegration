library(pbapply)
library(data.table)
library(rentrez)
set_entrez_key("49b3079321d573aaa12522e38a1b31d38e08")#ncbi account for dopreto to submit 10 queries per second
library(tidyverse)
library(ggplot2)

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
#[1] 20998

#search db
comention=pblapply(queries,function(x) {
		  request=entrez_search(db = "pubmed", term = x);
		  Sys.sleep(0.01);
		  return(request)})
#get succesful queries
cpgs=cbind(queries,sapply(comention,function(x) x[2]))
cpgs=cpgs[cpgs[,2]>0,]
cpgs=data.frame(cpgs[grep("AND  AND",cpgs[,1],invert=T),])
nrow(cpgs)
#[1] 707

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
#[1] 4018
tfs=data.frame(cbind(names(tfs),tfs))

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
mirs=data.frame(cbind(names(mirs),mirs))

#########################PER SUBTYPE#########################
#mix all together
pubmed=list(CpG=cpgs,TF=tfs,miRNA=mirs)

#BPs
length(unique(sapply(strsplit(unlist(lapply(pubmed,function(y) levels(y[,1]))),
							  " AND"),function(x) x[1])))
#[1] 146
#regulators
length(unique(sapply(strsplit(unlist(lapply(pubmed,function(y) levels(y[,1]))),
							  " AND"),function(x) x[2])))

#processes with the lots of regulators in PubMed have lots of genes?
inPM=table(sapply(strsplit(unlist(lapply(pubmed,function(y) levels(y[,1]))),
							  " AND"),function(x) x[1]))
entrez=unique(do.call(rbind,lapply(GS_GO_BP,function(x) 
	cbind(names(x),sapply(x,length)))))
entrez=entrez[entrez[,1]%in%names(inPM),]
entrez=entrez[order(match(entrez[,1],names(inPM))),]
entrez=data.frame(cbind(entrez,inPM),)
colnames(entrez)=c("process","genes","regulators")
entrez$genes=as.numeric(as.character(entrez$genes))
entrez$regulators=as.numeric(as.character(entrez$regulators))

png("inPMvsGenes.png")
 ggplot(genes,aes(x=genes,y=regulators,label=process))+geom_point()+
 ylab("regulators in PubMed")+theme(text=element_text(size=18))+
 geom_text(aes(label=ifelse(regulators>300,as.character(process),' '),
 	hjust=0.5,vjust=-0.3))
dev.off()

#recover queries per subtype
cpgs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="c"]))
cpgs=lapply(cpgs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	methy$UCSC_RefGene_Name[methy$IlmnID%in%x[[y]]],
	    	"(CpG methylation OR DNA methylation)",
	    	sep=" AND "))))
tfs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="E"]))
tfs=lapply(tfs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],
	    	myannot$hgnc_symbol[myannot$ensembl_gene_id%in%x[[y]]],
	    	sep=" AND "))))
mirs=lapply(regus,function(x) lapply(x,function(y) y[substr(y,1,1)=="h"]))
mirs=lapply(mirs,function(x) unlist(lapply(1:length(x),function(y) 
	    paste(names(x)[y],x[[y]],sep=" AND "))))
#mix all together
queries=lapply(1:5,function(x) c(cpgs[[x]],tfs[[x]],mirs[[x]]))
names(queries)=names(cpgs)
#how many queries per subtype were actually in pubmed?
queries=lapply(queries,function(x) sapply(pubmed,function(y) x[x%in%y[,1]]))
#data frame to plot
temp=data.frame(sapply(queries,function(x) sapply(x,length)))
temp=cbind(rownames(temp),temp)
colnames(temp)[1]="regulator"
temp=temp%>%pivot_longer(-regulator,names_to="subtype",values_to="sum")
temp$regulator=factor(temp$regulator,levels=c("CpG","TF","miRNA"))

png("InPubmed.png")
 ggplot(temp,aes(y=sum,x=subtype,fill=regulator))+
 geom_bar(position="dodge", stat="identity")+
 scale_fill_manual(values=gray.colors(5))+
 ylab("process AND regulator in pubmed")+xlab("")+
 theme(text=element_text(size=18))
dev.off()

#articles per omic
temp=data.frame(do.call(rbind,lapply(1:3,function(x) cbind(names(pubmed)[x],pubmed[[x]][,2]))))
colnames(temp)=c("regulator","articles")
temp$articles=as.numeric(as.character(temp$articles))
temp$regulator=factor(temp$regulator,levels=c("CpG","TF","miRNA"))

png("articles.png")
 ggplot(temp,aes(y=articles,x=regulator,fill=regulator))+geom_boxplot()+
 scale_fill_manual(values=gray.colors(5))+
 theme(text=element_text(size=18),legend.position="n")+xlab("")
dev.off()
