##############################################################################
library(biomaRt)
library(rentrez)
library(data.table)
library(pbapply)

##############################################################################
########### coherence with regulation models
##############################################################################
type=as.matrix(interacs[,3:4])
type[grep("^c",type[,1]),1]="CpG"
type[grep("^h",type[,1]),1]="miRNA"
type[type[,1]%in%TFtargets$TF,1]="TF"
type[!type[,1]%in%c("CpG","miRNA","TF"),1]="transcript"
potens=type[type$predictor!="transcript",]

#coeff distribution
png("coefficients")
 partR=ggplot(potens,aes(x=as.numeric(as.character(coef))))+geom_density(aes(group=predictor,color=predictor,fill=predictor),alpha=0.3)+xlim(-1,1)+scale_fill_manual(values=custom.col)+scale_color_manual(values=custom.col)+xlab("coefficient")
 partT=ggplot(type,aes(x=as.numeric(as.character(coef))))+geom_density(aes(group=predictor,color=predictor,fill=predictor),alpha=0.3)+xlim(-1,1)+scale_fill_manual(values=custom.col)+scale_color_manual(values=custom.col)+xlab("coefficient")
 grid.arrange(partR,partT,nrow=2)
dev.off()

png("coefficientsB.png")
 ggplot(type,aes(y=as.numeric(as.character(coef)),color=predictor))+geom_boxplot()+scale_color_manual(values=custom.col)+ylab("coefficient")
dev.off()

# CpGs enrichment for negative coeffs
length(intersect(grep("^c",interacs$predictor),which(interacs$coef<0)))
#[1] 5051
length(intersect(grep("^c",interacs$predictor),which(interacs$coef>0)))
#[1] 4532
length(intersect(grep("^c",interacs$predictor,invert=T),which(interacs$coef<0)))
#[1] 673
length(intersect(grep("^c",interacs$predictor,invert=T),which(interacs$coef>0)))
#[1] 1856
fisher.test(t(matrix(c(5051,673,4532,1856),ncol=2,nrow=2)),alternative="greater")
#p-value < 2.2e-16
# miRNAs enrichment for negative coeffs
length(intersect(grep("^h",interacs$predictor),which(interacs$coef<0)))
#[1] 33
length(intersect(grep("^h",interacs$predictor),which(interacs$coef>0)))
#[1] 76
length(intersect(grep("^h",interacs$predictor,invert=T),which(interacs$coef<0)))
#[1] 5691
length(intersect(grep("^h",interacs$predictor,invert=T),which(interacs$coef>0)))
#[1] 6312
fisher.test(t(matrix(c(33,5691,76,6312),ncol=2,nrow=2)),"greater")
#p-value = 0.0003342
# TFs enrichment for negative coeffs
sum(interacs$predictor%in%TFtargets$TF&interacs$coef<0)
#[1] 88
sum(interacs$predictor%in%TFtargets$TF&interacs$coef>0)
#[1] 224
length(intersect(intersect(which(!interacs$predictor%in%TFtargets$TF),grep("^h|c",interacs$predictor,invert=T,perl=T)),which(interacs$coef<0)))
#[1] 552
length(intersect(intersect(which(!interacs$predictor%in%TFtargets$TF),grep("^h|c",interacs$predictor,invert=T,perl=T)),which(interacs$coef>0)))
#[1] 1556
fisher.test(t(matrix(c(88,552,224,1556),ncol=2,nrow=2)))
#p-value = 0.4499

##############################################################################
########### am I linking genes with their KNOWN regulators?
##############################################################################
########### CpGs

#Input CpG-gene annotation = https://zwdzwd.github.io/InfiniumAnnotation
methy=fread("../ini/hm450.hg38.manifest.tsv",sep='\t',header=T,fill=T)
methy=methy[,c(5,21)]
methy=do.call(rbind,apply(methy,1,function(x) cbind(x[1],unlist(strsplit(x[2],";")))))

#build contingency tables for CpGs
#testMethyRegul=function(gen,subtype){
#anotadas=methy[which(methy[,2]==gen),1]
#seleccionadas=subtype$predictor[intersect(which(subtype$'pam50'==gen),grep("^c",subtype$predictor))]
#a=sum(anotadas%in%seleccionadas)
#b=sum(!anotadas%in%seleccionadas)
#c=sum(!seleccionadas%in%anotadas)
#d=length(unique(methy[,1]))-a-b-c
#signif=fisher.test(matrix(c(a,b,c,d),ncol=2,nrow=2),alternative="greater")$p.val
#I expect an overlap between annotated and selected regulators
#return(cbind(a,b,c,d,signif))}

#contingenciasMethy=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) testMethyRegul(y,x))))

#sign of the known CpGs 
nega=interacs[interacs$coef<0,]
negaCpG=nega[grep("^c",nega$predictor),]
negCpGanno=lapply(unique(negaCpG$predictor),function(x) sum(negaCpG$pam50[negaCpG$predictor==x]%in%methy[methy[,1]==x,2]))
table(unlist(negCpGanno))
#   0    1 
#4712    8 
posi=interacs[interacs$coef>0,]
posCpG=posi[grep("^c",posi$predictor),]
posCpGanno=lapply(unique(posCpG$predictor),function(x) sum(posCpG$pam50[posCpG$predictor==x]%in%methy[methy[,1]==x,2]))
table(unlist(posCpGanno))
#   0 
#4212 
fisher.test(matrix(c(8,0,4712,4212),ncol=2,nrow=2))
#p-value = 0.0085

#as a whole, there is no enrichment for CpGs
a=length(grep("^c",unique(interacs$predictor)))
c=384575-a
b=length(unique(interacs$predictor))-a
d=16475+433-b
CpGenri=t(matrix(c(a,c,b,d),ncol=2,nrow=2))
fisher.test(CpGenri,alternative="greater")$p.val
#p-value = 1

#there is neither an enrichment for known regulatory CpGs
cg=interacs[grep("^c",interacs$predictor),]
#a=sapply(unique(cg$predictor),function(x) sum(cg$pam50[cg$predictor==x]%in%methy[methy[,1]==x,2]))
#[1] 8
pam50=read.table("../ini/pam50.tsv",header=T)
methy[methy[,2]%in%pam50$hgnc_symbol,]
fisher.test(matrix(c(8,nrow(cg)-8,nrow(methy)-8,384575-nrow(cg)-nrow(methy)+8),ncol=2,nrow=2),alternative="greater")
#p-value = 1
fisher.test(matrix(c(8,9583-8,1182-8,384575-9583-1182+8),ncol=2,nrow=2),alternative="less")
#p-value = 2.447e-06

##############################################################################
########### miRNAs

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

#as a whole, there is an enrichment for miRs
a=length(grep("^h",unique(interacs$predictor)))
b=length(unique(interacs$predictor))-a
c=433-a
d=384575+16475-b
miRenri=t(matrix(c(a,c,b,d),ncol=2,nrow=2))
#      slctd NOslctd
#miR      97     336
#nomiR 10617  390433
fisher.test(miRenri,alternative="greater")
#p-value < 2.2e-16

##############################################################################
########### TFs

library(tftargets)#https://github.com/slowkow/tftargets

#transform TF target lists to tables easier to work with
load("TFtargets.RData")
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://jan2019.archive.ensembl.org")
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=95)
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

#only interested on TFs actually present on the input dataset
genes=rownames(DE.genes$her2_normal)
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name"), mart=mart)
myannot=myannot[myannot$ensembl_gene_id%in%genes,]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
TFtargets=TFtargets[TFtargets[,1]%in%myannot$hgnc_symbol,]

#test.TFregul=function(gen,subtype){
#anotadas=unique(TFtargets[which(TFtargets$'target'==gen),1])
#seleccionadas=as.character(subtype$predictor[intersect(which(subtype$'pam50'==gen),which(subtype$predictor%in%TFtargets$TF))])
#a=sum(anotadas%in%seleccionadas)
#b=sum(!anotadas%in%seleccionadas)
#c=sum(!seleccionadas%in%anotadas)
#d=length(unique(TFtargets$TF))-a-b-c
#signif=fisher.test(matrix(c(a,b,c,d),ncol=2,nrow=2),alternative="greater")$p.val
#return(cbind(a,b,c,d,signif))}

#contingencias.ENSG=pblapply(interacs,function(x) t(sapply(unique(x$pam50),function(y) test.TFregul(y,x))))
#genes with ENSGs selected
#validatedTFs=do.call(rbind,lapply(c(1,3:5),function(x) 
#	cbind(names(interacs)[x],as.character(unique(interacs[[x]]$pam50)[contingencias.ENSG[[x]][,1]>0]),contingencias.ENSG[[x]][contingencias.ENSG[[x]][,1]>0,])))

#as a whole, there is an enrichment for TFs
a=sum(unique(interacs$predictor)%in%TFtargets$TF)
b=length(unique(interacs$predictor))-a
c=length(unique(TFtargets$TF))-a
d=384575+433+16475-length(unique(interacs$predictor))
TFenri=t(matrix(c(a,c,b,d),ncol=2,nrow=2))
#     slctd NOslctd
#TF     268    2011
#noTF 10446  390769
fisher.test(TFenri,alternative="greater")
#p-value < 2.2e-16

#there is also an enrichment for known TFs
tfs=interacs[interacs$predictor%in%TFtargets$TF,]
tfs=lapply(unique(tfs$subtype),function(x) tfs[tfs$subtype==x,])
a=lapply(tfs,function(y) sapply(unique(y$predictor),function(x) 
	sum(y$pam50[y$predictor==x]%in%interesTF$target[interesTF$TF==x])))
sum(sapply(a,sum))
#[1] 41
interesTF=TFtargets[TFtargets$target%in%pam50$hgnc_symbol,]
fisher.test(matrix(c(41,sum(sapply(tfs,nrow))-41,nrow(interesTF)-41,nrow(TFtargets)-sum(sapply(tfs,nrow))-nrow(interesTF)+41),ncol=2,nrow=2),alternative="greater")
#p-value < 2.2e-16

#####################################################################
#######how many times terms are comentioned in literature
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
######coherence of regulation
#####################################################################
#negative coeff CpGs methylation is opposite to target gene expression	     
cpgs=interacs[grep("^c",interacs$predictor),]
cpgs=lapply(unique(cpgs$subtype),function(x) cpgs[cpgs$subtype==x,])
dvalues=lapply(1:4,function(x) t(apply(cpgs[[x]],1,function(z) 
	rbind(DE.genes[[x]]$logFC[rownames(DE.genes[[x]])==z[2]],DM.cpg1[[x]]$logFC[rownames(DM.cpg1[[x]])==z[3]]))))
sameSign=lapply(dvalues,function(x) rowSums(x<0))
sameSign=lapply(1:4,function(x) cbind(cpgs[[x]]$coef,sameSign[[x]]))
table(unlist(lapply(1:4,function(x) sameSign[[x]][sameSign[[x]][,1]<0,2])))
#   0    1    2 
#1083 1864  637
table(unlist(lapply(1:4,function(x) sameSign[[x]][sameSign[[x]][,1]>0,2])))
#   0    1    2 
#1002 2089 1147
fisher.test(matrix(c(1864,2089,1720,2149),ncol=2,nrow=2),alternative="less")
#p-value = 0.008846

#negative coeff miRs expression is NOT opposite to target gene expression	     
mirs=interacs[grep("^h",interacs$predictor),]
mirs=lapply(unique(mirs$subtype),function(x) mirs[mirs$subtype==x,])
dvalues.miR=lapply(1:4,function(x) t(apply(mirs[[x]],1,function(z) 
	rbind(DE.genes[[x]]$logFC[rownames(DE.genes[[x]])==z[2]],DE.miR[[x]]$logFC[rownames(DE.miR[[x]])==z[3]]))))
sameSign.miR=lapply(dvalues.miR,function(x) rowSums(x<0))
sameSign.miR=lapply(1:4,function(x) cbind(mirs[[x]]$coef,sameSign.miR[[x]]))
table(unlist(lapply(1:4,function(x) sameSign.miR[[x]][sameSign.miR[[x]][,1]<0,2])))
#0 1 2 
#9 8 1 
table(unlist(lapply(1:4,function(x) sameSign.miR[[x]][sameSign.miR[[x]][,1]>0,2])))
# 0  1  2 
#21 24  6 
fisher.test(matrix(c(8,24,10,27),ncol=2,nrow=2))
#p-value = 1

#negative coeff TFs expression is NOT opposite to target gene expression	     
tfs=interacs[interacs$predictor%in%TFtargets$TF,]
tfs=lapply(unique(tfs$subtype),function(x) tfs[tfs$subtype==x,])
dvalues.TF=lapply(1:4,function(x) t(apply(tfs[[x]],1,function(z) 
	rbind(DE.genes[[x]]$logFC[rownames(DE.genes[[x]])==z[2]],DE.genes[[x]]$logFC[rownames(DE.genes[[x]])==z[3]]))))
sameSign.TF=lapply(dvalues.TF,function(x) rowSums(x<0))
sameSign.TF=lapply(1:4,function(x) cbind(tfs[[x]]$coef,sameSign.TF[[x]]))
table(unlist(lapply(1:4,function(x) sameSign.TF[[x]][sameSign.TF[[x]][,1]<0,2])))
# 0  1  2 
#19 30 21 
table(unlist(lapply(1:4,function(x) sameSign.TF[[x]][sameSign.TF[[x]][,1]>0,2])))
#0  1  2 
#54 86 53 
fisher.test(matrix(c(30,86,40,107),ncol=2,nrow=2))
#p-value = 0.8884

#joining all regulators negative coeff predictors are NOT opposite to target gene expression	     
colSums(rbind(table(unlist(lapply(1:4,function(x) sameSign[[x]][sameSign[[x]][,1]<0,2]))),table(unlist(lapply(1:4,function(x) sameSign.miR[[x]][sameSign.miR[[x]][,1]<0,2]))),table(unlist(lapply(1:4,function(x) sameSign.TF[[x]][sameSign.TF[[x]][,1]<0,2])))))
#   0    1    2 
#1030 2127 1169 
colSums(rbind(table(unlist(lapply(1:4,function(x) sameSign[[x]][sameSign[[x]][,1]>0,2]))),table(unlist(lapply(1:4,function(x) sameSign.miR[[x]][sameSign.miR[[x]][,1]>0,2]))),table(unlist(lapply(1:4,function(x) sameSign.TF[[x]][sameSign.TF[[x]][,1]>0,2])))))
#   0    1    2 
#1158 1974  696 
fisher.test(matrix(c(2127,1974,2199,1854),ncol=2,nrow=2),alternative="less")
#p-value = 0.01615