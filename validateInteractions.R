##############################################################################
########### methy: am I linking genes with their KNOWN regulatory CpGs?
##############################################################################
#Input CpG-gene annotation = https://zwdzwd.github.io/InfiniumAnnotation
methy=fread("ini/hm450.hg38.manifest.tsv",sep='\t',header=T,fill=T)
methy=methy[,c(5,21)]
methy=do.call(rbind,apply(methy,1,function(x) cbind(x[1],unlist(strsplit(x[2],";")))))

#build contingency tables for CpGs
testMethyRegul=function(gen,subtype){
anotadas=methy[which(methy[,2]==gen),1]
seleccionadas=subtype$predictor[intersect(which(subtype$'pam50'==gen),grep("^c",subtype$predictor))]
a=sum(anotadas%in%seleccionadas)
b=sum(!anotadas%in%seleccionadas)
c=sum(!seleccionadas%in%anotadas)
d=length(unique(methy[,1]))-a-b-c
signif=fisher.test(matrix(c(a,b,c,d),ncol=2,nrow=2),alternative="greater")$p.val
#I expect an overlap between annotated and selected regulators
return(cbind(a,b,c,d,signif))}

contingenciasMethy=pblapply(interacs,function(x) sapply(unique(x$pam50),function(y) t(testMethyRegul(y,x))))
#genes with CpGs selected
sapply(contingencias,function(x) sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#       44        42        45        44        45 
sapply(contingenciasMethy,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#proportion of pam50 genes per subtype with known CpG regulators selected
#     Basal       Her2       LumA       LumB  non-tumor 
#0.02272727 0.00000000 0.08888889 0.00000000 0.00000000 
#sampe proportions when I test alternative="both"
#sapply(contingenciasMethy,function(x) sum(apply(x[,1:4],1,function(y) fisher.test(matrix(y,ncol=2,nrow=2))$p.val)<0.01)/nrow(x))     

##############################################################################
########### miR: am I linking genes with their KNOWN regulatory miRr?
##############################################################################
library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html
#Searching mirecords, mirtarbase, tarbase,diana_microt,elmmo, microcosm, miranda, mirdbpictar, pita, targetscan, pharmaco_mir ...

#get all validated miRs interactions for PAM50 genes
miRtargets=get_multimir(target=pam50$ensembl_gene_id,summary=F,table="all",legacy.out=F)
miRtargetsV=select(miRtargets,keys="validated",columns=columns(miRtargets),keytype="type")
miRtargetsP=select(miRtargets,keys="predicted",columns=columns(miRtargets),keytype="type")
#i want miR stem names coz my datset miR names are for stem stage  
miRtargetsV=cbind(miRtargetsV,sapply(strsplit(miRtargetsV$mature_mirna_id,"-"),function(x) paste(x[1],x[2],x[3],sep='-')))

#I'm only interested on the interacs I could find with my miR dataset
load("../ini/porSubti.RData")
miR=rownames(concatenadas$Her2$mir)
miRtargetsV=miRtargetsV[miRtargetsV[,14]%in%miR|miRtargetsV$mature_mirna_id%in%miR,]

test.miRegul=function(gen,subtype,annot){
anotadas=unique(unlist(annot[which(annot$target_symbol==gen),c(3,14)]))
seleccionadas=subtype$predictor[intersect(which(subtype$'pam50'==gen),grep("^h",subtype$predictor))]
a=sum(anotadas%in%seleccionadas)
b=sum(!anotadas%in%seleccionadas)
c=sum(!seleccionadas%in%anotadas)
d=length(miR)-a-b-c
signif=fisher.test(matrix(c(a,b,c,d),ncol=2,nrow=2),alternative="greater")$p.val
return(cbind(a,b,c,d,signif))}
contingencias.miR=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) test.miRegul(y,x,miRtargetsV))))
#genes with miRs selected
sapply(contingencias.miR,function(x) sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#        9         3        12        13         8 
#models didnt selected regulataroy miR
sapply(contingencias.miR,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#        0         0         0         0         0 

#predicted miR are neither selected
miRtargetsP=cbind(miRtargetsP,sapply(strsplit(miRtargetsP$mature_mirna_id,"-"),function(x) paste(x[1],x[2],x[3],sep='-')))
miRtargetsP=miRtargetsP[miRtargetsP[,14]%in%miR|miRtargetsP$mature_mirna_id%in%miR,]
contingencias.miR=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) test.miRegul(y,x,miRtargetsP))))
sapply(contingencias.miR,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#        0         0         0         0         0 

#but there are very few interactions that I could have found
dim(miRtargetsP)
#[1] 89 14
dim(miRtargetsV)
#[1] 32 14

##############################################################################
########### TFs: am I linking genes with their KNOWN tfs?
##############################################################################
##############################################################################
########### TFs
##############################################################################
library(tftargets)#https://github.com/slowkow/tftargets

#transform TF target lists to tables easier to work with


TF=read.table("data/TFCheckpoint_download_180515.txt",header=T,sep='\t',fill=T,quote="")
temp=myannot[myannot$ensembl_gene_id%in%unique(interacs[,2]),c(1,7)]
i=which(interacs[,2]%in%temp$ensembl_gene_id[temp$entrezgene%in%TF$entrez_human])
length(i)
#[1] 3553 interactions with TFs on TFCheckpoint

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

intrcsTF=apply(interacs[i,4:5],1,function(x) knownTF(x[2],x[1]))
sum(intrcsTF!="")
#[1] 202 interactions with TFs reported in DB
interacs=cbind(interacs,NA)
interacs[i,6]=intrcsTF
colnames(interacs)=c("pam50","predictor","coef","pam50Symbol","predictorSymbol","TFsupportedBy")
length(unique(interacs[interacs[,6]!=""&!is.na(interacs[,6]),5]))
#[1] 131 TFs with known TF-target interaction in tftargets



lambda=unlist(lapply(quality,function(x) as.numeric(x[,3])))
lambda=as.data.frame(cbind(lambda,as.character(sapply(names(quality),rep,50))))
colnames(lambda)[2]="subtype"
ggplot(lambda,aes(x=as.numeric(lambda)))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+scale_x_continuous(minor_breaks=NULL)
rmse=unlist(lapply(quality,function(x) as.numeric(x[,4])))
rmse=as.data.frame(cbind(rmse,as.character(sapply(names(quality),rep,50))))
colnames(rmse)[2]="subtype"
ggplot(rmse,aes(x=as.numeric(rmse)))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+scale_x_continuous(minor_breaks=NULL)

library(biomaRt)
library(rentrez)

mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://apr2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","external_gene_name","entrezgene"), mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]

#hgnc_symbol instead of ensemblID
interacs=unique(do.call(rbind,sifs))
interacs=interacs[order(interacs[,1]),]
interacs=interacs[interacs[,2]!="(Intercept)",]
targets=table(interacs[,1])
temp=sapply(1:50,function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(targets)[x]],targets[x]))
interacs=cbind(interacs,unlist(temp))
interacs=interacs[order(interacs[,2]),]
i=grep("ENSG",interacs[,2])
interacs=cbind(interacs,interacs[,2])
notarg=table(interacs[i,2])
temp=sapply(1:length(notarg),function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(notarg)[x]],notarg[x]))
interacs[i,5]=unlist(temp)

#####################################################################
#######how many times terms are mentioned in literature
#####################################################################
#both on the same paper
comention=pbsapply(round(seq(1,11119,length=20)),function(i)
  apply(interacs[i:(i+585),4:5],1,function(x) {
  reque=entrez_search(db = "pubmed", term = paste(x[1]," AND ",x[2],collapse=" "));
  Sys.sleep(0.1);
 return(reque)}))

cuentas=as.numeric(apply(comention,c(1,2),function(x) unlist(x[[1]][2])))
query=as.character(apply(comention,c(1,2),function(x) unlist(x[[1]][4])))
comen=cbind(query,cuentas)
comen=comen[comen[,2]!="0",]
#cat temp|perl -pe 'unless(/ AND \(/){s/.*?\sOR//g}unless(/\) AND \(/){s/\(.*?OR//g;s/\".*?OR //;s/(^ |\"|\)|\[All Fields\])//g}'
comention=read.table("Downloads/temp",sep='\t',header=T)
temp=t(sapply(strsplit(as.character(comention$query)," +AND +",perl=T),as.character))
comention=cbind(comention,temp)
comention[,3]=toupper(comention[,3])
comention[,4]=toupper(comention[,4])

#each mentioned
query=unique(unlist(comention[,3:4]))
mention=pbsapply(query,function(x){
	reque=entrez_search(db = "pubmed", term = x);
	Sys.sleep(0.1);
	return(reque)})
mention=as.matrix(t(mention)[,2])

i=apply(comention,1,function(x) which(interacs[,3]==x[3]&interacs[,4]==x[4]))
interacs=cbind(interacs,NA)
interacs[i,7]=comention$cuentas
temp=t(apply(interacs[i,],1,function(x) 
	cbind(mention[rownames(mention)==x[3],2],mention[rownames(mention)==x[4],2])))
interacs=cbind(interacs,NA,NA)
interacs[i,8]=temp[,1]
interacs[i,9]=temp[,2]
colnames(interacs)[7:9]=c("comention","pam50Mention","predictorMention")

#####################################################################
######all together per subtype
#####################################################################

sifs1=lapply(sifs,function(x) x[x[,2]!="(Intercept)",])
temp=lapply(sifs,function(x) do.call(rbind,apply(x,1,function(y) interacs[interacs[,1]==y[1]&interacs[,2]==y[2],])))
sifs1=lapply(1:5,function(x) cbind(as.matrix(sifs1[[x]][,3]),as.matrix(temp[[x]])))
names(sifs1)=names(sifs)
sifs1=lapply(sifs1,function(x) x[,c(2:3,1,4:ncol(x))])
#####################################################################
######coherence of regulation
#####################################################################
	     
tempi=sifs1[c(1:2,5,3)]
i.cpg=lapply(tempi,function(x) grep("hsa|ENSG",x[,2],perl=T,invert=T))
i.mir=lapply(tempi,function(x) grep("hsa",x[,2]))

DA=lapply(1:4,function(x) rbind(DE.genes[[x]],DE.miR[[x]],DM.cpg1[[x]]))
up=lapply(1:4,function(x) sapply(tempi[[x]][,1],function(y) DE.genes[[x]]$logFC[rownames(DE.genes[[x]])==y]))
prediUP=lapply(1:4,function(x) sapply(tempi[[x]][,2],function(y) DA[[x]]$logFC[rownames(DA[[x]])==y]))
tempi=lapply(1:4,function(x) cbind(tempi[[x]],up[[x]]>0,prediUP[[x]]>0))
#lapply(1:4,function(x) fisher.test(table(as.data.frame(tempi[[x]][i.mir[[x]],8:9]))))#no significant
tempi$normal=sifs1$normal
sifs1=tempi
colnames(sifs1[[1]])[c(3,10:11)]=c("beta","pam50UP","predictorUP")
colnames(sifs1[[2]])[c(3,10:11)]=c("beta","pam50UP","predictorUP")
colnames(sifs1[[3]])[c(3,10:11)]=c("beta","pam50UP","predictorUP")
colnames(sifs1[[4]])[c(3,10:11)]=c("beta","pam50UP","predictorUP")
colnames(sifs1[[5]])[3]="beta"
save(sifs1,file="annotatedSifsAlpha0.5.RData")
