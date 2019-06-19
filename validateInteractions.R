library(biomaRt)
library(rentrez)
library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name","entrezgene"), mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]

#hgnc_symbol instead of ensemblID
interacs=unique(do.call(rbind,sifs))
interacs=interacs[order(interacs[,1]),]
interacs=interacs[interacs[,2]!="",]
interacs=interacs[interacs[,2]!="(Intercept)",]
targets=table(interacs[,1])
temp=sapply(1:50,function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(targets)[x]],targets[x]))
interacs=cbind(interacs,unlist(temp))
interacs=interacs[order(interacs[,2]),]
i=grep("ENSG",interacs[,2])
interacs=cbind(interacs,interacs[,2])
notarg=table(interacs[i,4])
temp=sapply(1:length(notarg),function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(notarg)[x]],notarg[x]))
interacs[i,4]=unlist(temp)

#how many times both entities are mentioned in the same paper
comention=apply(interacs[,3:4],1,function(x)  
	rentrez::entrez_search(db = "pubmed", term = paste(x,collapse=" "))$count)

##############################################################################
########### methy
##############################################################################

methy=read.table("ini/hm450.hg38.manifest.tsv",sep='\t',header=T,fill=T)
#am I linking genes with their known regulating CpGs? NOP
i=grep("ENSG|hsa",interacs[,2],perl=T,invert=T)
length(i)
#[1] 3997
knownReg=sapply(i,function(x) 
	grep(interacs[x,3],methy$gene_HGNC[methy$probeID==interacs[x,4]]))
sum(sapply(knownReg,length)==0)
#[1] 3997
#are they in the same chr?
sameChr=apply(interacs[i,1:2],1,function(x) 
	sum(myannot$chromosome_name[myannot$ensembl_gene_id==x[1]]==methy$CpG_chrm[methy$probeID==x[2]]))
sum(sameChr!=0) 
#[1] 156 sólo estos están en el mismo cromosoma y están lejos ↓
interacs=cbind(interacs,NA)
interacs[i,5]="n"
interacs[i[sameChr!=0],5]="y"
temp=interacs[interacs[,5]=="y",1:2]#todos afectan PTTG1  
apply(methy[methy$probeID%in%temp[,2],2:3]-myannot$start_position[myannot$ensembl_gene_id=="ENSG00000164611"],2,function(x) min(abs(x)))
#CpG_beg CpG_end 
#4076136 4076138 
lapply(sifs,function(x) sum(x[,1]=="ENSG00000164611"&x[,2]%in%temp[,2]))
$Basal
[1] 54
$Her2
[1] 72
$LumB
[1] 30
$normal
[1] 0

##############################################################################
########### miR
##############################################################################
#get all the interactions of selected miRs
mirInteractions=do.call(rbind,lapply(unique(as.character(interacs[i,2])),function(x) 
	get_multimir(mirna = x, summary = F)@data))

#found interactions have been reported? NOP
i=grep("hsa",interacs[,2])
knownTarget=apply(interacs[i,],1,function(x) 
	sum(mirInteractions$mature_mirna_id==x[2]&mirInteractions$target_symbol==x[3]))
sum(knownTarget)
#[1] 0

##############################################################################
########### TFs
##############################################################################
