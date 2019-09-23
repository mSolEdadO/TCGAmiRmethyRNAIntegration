##############################################################################
library(biomaRt)
library(rentrez)
library(data.table)
library(pbapply)

########### methy: am I linking genes with their KNOWN regulatory CpGs?
##############################################################################
#Input CpG-gene annotation = https://zwdzwd.github.io/InfiniumAnnotation
methy=fread("../ini/hm450.hg38.manifest.tsv",sep='\t',header=T,fill=T)
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

contingenciasMethy=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) testMethyRegul(y,x))))
#genes with CpGs selected
sapply(contingenciasMethy,function(x) sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#       44        42        45        44        45
#known CpG interactios from the total selected per subtype 
sapply(contingenciasMethy,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#       Basal         Her2         LumA         LumB    non-tumor 
#0.0004214075 0.0000000000 0.0019002375 0.0000000000 0.0000000000 
#sampe proportions when I test alternative="both"
#sapply(contingenciasMethy,function(x) sum(apply(x[,1:4],1,function(y) fisher.test(matrix(y,ncol=2,nrow=2))$p.val)<0.01)/nrow(x))     

##############################################################################
########### miR: am I linking genes with their KNOWN regulatory miRr?
##############################################################################
library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html
#Searching mirecords, mirtarbase, tarbase,diana_microt,elmmo, microcosm, miranda, mirdbpictar, pita, targetscan, pharmaco_mir ...
pam50=read.table("../ini/pam50.tsv",header=T)

#get all validated miRs interactions for PAM50 genes
miRtargets=get_multimir(target=pam50$ensembl_gene_id,summary=F,table="all",legacy.out=F)
miRtargetsV=select(miRtargets,keys="validated",columns=columns(miRtargets),keytype="type")
miRtargetsP=select(miRtargets,keys="predicted",columns=columns(miRtargets),keytype="type")
#i want miR stem names coz my datset miR names are for stem stage  
miRtargetsV=cbind(miRtargetsV,sapply(strsplit(miRtargetsV$mature_mirna_id,"-"),function(x) 
	paste(x[1],x[2],x[3],sep='-')))

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
#interactions with miR selected per subtype
sapply(contingencias.miR,function(x) sum(x[,3]))
#    Basal      Her2      LumA      LumB non-tumor 
#       15         4        20        30        40 
#models didnt selected regulataroy miR
sapply(contingencias.miR,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#        0         0         0         0         0 

#predicted miR are neither selected
miRtargetsP=cbind(miRtargetsP,sapply(strsplit(miRtargetsP$mature_mirna_id,"-"),function(x)
 paste(x[1],x[2],x[3],sep='-')))
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
library(tftargets)#https://github.com/slowkow/tftargets

#transform TF target lists to tables easier to work with
load("TFtargets.RData")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://jan2019.archive.ensembl.org")
myannot=getBM(
 attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","external_gene_name","entrezgene"), 
 mart=mart)
myannot=myannot[!is.na(myannot$entrezgene),]
myannot=myannot[!duplicated(myannot$entrezgene),]
#get hgnc_id for data annotated with entrez id
TFtargets$TRED=TFtargets$TRED[TFtargets$TRED[,2]%in%myannot$entrezgene,]
TFtargets$TRED=TFtargets$TRED[order(TFtargets$TRED[,2]),]
temp=table(TFtargets$TRED[,2])
ids=unlist(sapply(1:length(temp),function(x)
 rep(as.character(myannot$hgnc_symbol)[which(as.character(myannot$entrezgene)==names(temp)[x])],temp[x])))
TFtargets$TRED[,2]=unlist(ids)

#only interested on TFs actually presente on the input dataset
genes=rownames(DE.genes$her2_normal)
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name"), mart=mart)
myannot=myannot[myannot$ensembl_gene_id%in%genes,]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
TFtargets=TFtargets[TFtargets[,1]%in%myannot$hgnc_symbol,]

test.TFregul=function(gen,subtype){
anotadas=unique(TFtargets[which(TFtargets$'target'==gen),1])
seleccionadas=as.character(subtype$predictor[intersect(which(subtype$'pam50'==gen),which(subtype$predictor%in%TFtargets$TF))])
a=sum(anotadas%in%seleccionadas)
b=sum(!anotadas%in%seleccionadas)
c=sum(!seleccionadas%in%anotadas)
d=length(unique(TFtargets$TF))-a-b-c
signif=fisher.test(matrix(c(a,b,c,d),ncol=2,nrow=2),alternative="greater")$p.val
return(cbind(a,b,c,d,signif))}

contingencias.ENSG=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) test.TFregul(y,x))))
#genes with ENSGs selected
sapply(contingencias.ENSG,function(x) sum(x[,3]>0))
#    Basal      Her2      LumA      LumB non-tumor 
#       36        19        40        32        39 
#interactions with ENSG selected per subtype
#sapply(contingencias.ENSG,function(x) sum(x[,3]))
#    Basal      Her2      LumA      LumB non-tumor 
#      369        73      1138       341       458 
#interactions with TFs selected per subtype
sapply(contingencias.ENSG,function(x) sum(x[,3]))
#    Basal      Her2      LumA      LumB non-tumor 
#       40         9       140        40        42 
#predicted TFs are neither selected
sapply(contingencias.ENSG,function(x) sum(x[,5]<0.01)/sum(x[,3]>0))
#     Basal       Her2       LumA       LumB  non-tumor 
#0.00000000 0.00000000 0.05000000 0.00000000 0.02564103 
validatedTFs=do.call(rbind,lapply(c(1,3:5),function(x) 
	cbind(names(interacs)[x],as.character(unique(interacs[[x]]$pam50)[contingencias.ENSG[[x]][,1]>0]),contingencias.ENSG[[x]][contingencias.ENSG[[x]][,1]>0,])))

#as a whole, there is no enrichment for TFs
coefs=fread("slctdPrdctrs.csv",header=T)
coefs=coefs[grep("^h|c",as.character(coefs$predictor),perl=T,invert=T),]
sum(unique(TFtargets$TF)%in%unique(coefs$predictor))
#[1] 274
sum(!unique(TFtargets$TF)%in%unique(coefs$predictor))
#[1] 2041
sum(!unique(coefs$predictor)%in%unique(TFtargets$TF))
#[1] 1858
fisher.test(matrix(c(274,2041,1858,16475-274-2041-1858),nrow=2,ncol=2))$p.val
#[1] 0.9603

#####################################################################
#######how many times terms are mentioned in literature
#####################################################################
#both on the same paper
set_entrez_key("49b3079321d573aaa12522e38a1b31d38e08")#ncbi account for dopreto to submit 10 queries per second

comention=pbsapply(seq(1,12112,10),function(i) 
	apply(interacs[i:(i+9),],1,function(x) {
		reque=entrez_search(db = "pubmed", term = paste(x[2]," AND ",x[3],collapse=" "));
		Sys.sleep(0.1);
		return(reque)}))


cuentas=as.numeric(apply(comention,c(1,2),function(x) unlist(x[[1]][2])))
comen=cbind(interacs[,2:3],cuentas[1:12112])#further than interacs nrow, search was NA AND NA
comen=comen[comen[,2]!="0",]
colnames(comen)[3]="comention"
write.table(comen,"comention.tsv",sep='\t',quote=F,row.names=F)

#each mentioned
query=unique(unlist(comen[,3:4]))
mention=pbsapply(seq(1,10746,7),function(i)
	sapply(query[i:(i+6)],function(x){
	reque=entrez_search(db = "pubmed", term = x);
	Sys.sleep(0.1);
	return(reque)}))
#gives a matrix of output lines(5)*searched terms(7) rows per length(seq(1,10746,7)) columns
cuentas=unlist(apply(mention,2,function(x) x[seq(2,35,5)]))[1:length(query)]#2 is the index of id counts in the columns
#further than query len, search was NA
mention=cbind(query,cuentas)

continComent=function(interac){
	members=as.character(interac)
	a=interac$comention
	b=mention$cuentas[mention$query==members[1]]-a
	c=mention$cuentas[mention$query==members[2]]-a
	d=30151833-b-c-a
	mat=matrix(c(a,b,c,d),ncol=2,nrow=2)
	pval=fisher.test(mat,alternative="greater")$p.val
	return(c(mat,pval))}
fshr=sapply(1:nrow(comen),function(x) continComent(comen))
fshr=cbind(fshr,p.adjust(fshr[,5],"fdr")
colnames(fshr)=c("comention","pam50mention","predimention","neither","p.val","fdr")
comention=cbind(comention,fshr)
write.table(comention,"comention.tsv",sep='\t',quote=F,row.names=F)

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

lambda=unlist(lapply(quality,function(x) as.numeric(x[,3])))
lambda=as.data.frame(cbind(lambda,as.character(sapply(names(quality),rep,50))))
colnames(lambda)[2]="subtype"
ggplot(lambda,aes(x=as.numeric(lambda)))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+scale_x_continuous(minor_breaks=NULL)
rmse=unlist(lapply(quality,function(x) as.numeric(x[,4])))
rmse=as.data.frame(cbind(rmse,as.character(sapply(names(quality),rep,50))))
colnames(rmse)[2]="subtype"
ggplot(rmse,aes(x=as.numeric(rmse)))+geom_density(aes(group=subtype,color=subtype,fill=subtype),alpha=0.3)+scale_x_continuous(minor_breaks=NULL)
