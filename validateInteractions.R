library(biomaRt)
library(rentrez)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","external_gene_name","entrezgene"), mart=mart)
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
#[1] 4202
knownReg=sapply(i,function(x) 
	grep(interacs[x,3],methy$gene_HGNC[methy$probeID==interacs[x,4]]))
sum(sapply(knownReg,length)==0)
#[1] 4202
#are they in the same chr?
sameChr=apply(interacs[i,1:2],1,function(x) 
	sum(myannot$chromosome_name[myannot$ensembl_gene_id==x[1]]==methy$CpG_chrm[methy$probeID==x[2]]))
sum(sameChr!=0) 
#[1] 173 sólo estos están en el mismo cromosoma y están lejos ↓
interacs=cbind(interacs,NA)
interacs[i,6]="n"
interacs[i[sameChr!=0],6]="y"
temp=interacs[which(interacs[,6]=="y"),]#todos afectan PTTG1  
apply(methy[methy$probeID%in%temp[,2],2:3]-myannot$start_position[myannot$ensembl_gene_id=="ENSG00000164611"],2,function(x) min(abs(x)))
#CpG_beg CpG_end 
#2846047 2846045 
lapply(sifs,function(x) sum(x[,1]=="ENSG00000164611"&x[,2]%in%temp[,2]))
$Basal
[1] 54
$Her2
[1] 72
$LumB
[1] 30
$normal
[1] 0
$LumA
[1] 17

##############################################################################
########### miR
##############################################################################
library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html
#Searching mirecords, mirtarbase, tarbase,diana_microt,elmmo, microcosm, miranda, mirdbpictar, pita, targetscan, pharmaco_mir ...

#get all the interactions of selected miRs
i=grep("hsa",interacs[,2])
length(unique(interacs[i,2]))
#[1] 241
miRNAs   = list_multimir("mirna")
sum(unique(interacs[i,2])%in%miRNAs$mature_mirna_id)
#[1] 7
mirInteractions=lapply(miRNAs$mature_mirna_id[miRNAs$mature_mirna_id%in%unique(interacs[i,2])],function(x)
	get_multimir(mirna = x, summary = F,table="all")@data)
#found interactions have been reported? NOP
knownTarget=sapply(mirInteractions,function(x)
	sum(interacs[i[interacs[i,2]%in%x$mature_mirna_id],3]%in%x$target_symbol))
sum(knownTarget)
#[1] 0
mirInteractions=lapply(unique(interacs[i,3]),function(x)
	get_multimir(target = x, summary = F,table="all")@data)
knownTarget=sapply(mirInteractions,function(x)
	sum(interacs[i[interacs[i,3]%in%x$target_symbol],2]%in%x$mature_mirna_id))
sum(knownTarget)
#[1] 0

##############################################################################
########### TFs
##############################################################################
library(tftargets)#https://github.com/slowkow/tftargets

#TF=read.table("data/TFCheckpoint_download_180515.txt",header=T,sep='\t',fill=T,quote="")
#temp=myannot[myannot$ensembl_gene_id%in%unique(interacs[,2]),c(1,4)]
#i=which(interacs[,2]%in%temp$ensembl_gene_id[temp$entrezgene%in%TF$entrez_human])
#length(i)
#[1] 2135 interactions with TFs

pam50annot=myannot[myannot$ensembl_gene_id%in%unique(interacs[,1]),]
i=grep("ENSG",interacs[,2])

knownTF=function(TF,target){
support=character()
if(pam50annot$entrezgene[pam50annot$hgnc_symbol==target]%in%TRED[[TF]]){
	support=c(support,"TRED")}#Predicted and known targets
if(target%in%ITFP[[TF]]){support=c(support,"ITFP")} #predicted
if(pam50annot$entrezgene[pam50annot$hgnc_symbol==target]%in%ENCODE[[TF]]){
	support=c(support,"ENCODE")}#ChipSeq
if(target%in%Neph2012[[TF]]){support=c(support,"Neph2012")}#DNAseq footprinting
if(target%in%TRRUST[[TF]]){support=c(support,"TRRUST")}#small-scale experiments
if(target%in%Marbach2016[[TF]]){support=c(support,"Marbach2016")}
return(paste(support,collapse=','))}#CAGE+binding motifs
#no sé cuando recopilaron los datos

intrcsTF=apply(interacs[i,3:4],1,function(x) knownTF(x[2],x[1]))
sum(intrcsTF!="")
#[1] 188 interactions with TFs reported in DB
interacs=cbind(interacs,NA)
intercs[i,6]=intrcsTF
colnames(interacs)=c("pam50","predictor","pam50Symbol","predictorSymbol","TFsupportedBy")
length(unique(interacs[interacs[,5]!=""&!is.na(interacs[,5]),4]))
#[1] 131 TFs involved in known TF-target interaction
