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
########### CpGs

#Input CpG-gene annotation = https://zwdzwd.github.io/InfiniumAnnotationZ<<<<<<<<<--NO
methy=fread("HumanMethylation450_15017482_v1-2.csv",sep=',',header=T,fill=T)#use Illumina file
methy=methy[,c(1,22,24)]
methy=methy[methy[,2]!="",]
methy=do.call(rbind,apply(methy,1,function(x) 
	cbind(x[1],unlist(strsplit(x[2],";")),unlist(strsplit(x[3],";")))))
writr.table(methy,"MapMethy.tsv",sep='\t',quote=F,row.names=F)

withCpG=do.call(rbind,lapply(unique(interacs$pam50Symbol),function(x)
	interacs[interacs$pam50Symbol==x,c(1,3:5)][which(interacs$predictor[interacs$pam50Symbol==x]%in%methy[methy[,2]==x,1]),]))

##############################################################################
########### miRNAs

library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html
#Searching mirecords, mirtarbase, tarbase,diana_microt,elmmo, microcosm, miranda, mirdbpictar, pita, targetscan, pharmaco_mir ...
pam50=read.table("../ini/pam50.tsv",header=T)

#get all validated miRs interactions for PAM50 genes
miRtargets=get_multimir(target=pam50$ensembl_gene_id,summary=F,table="all",legacy.out=F)
miRtargetsV=select(miRtargets,keys="validated",columns=columns(miRtargets),keytype="type")
miRtargetsP=select(miRtargets,keys="predicted",columns=columns(miRtargets),keytype="type")
mirIDs=fread(miR.ids.map.tsv)#from mirbase <---------LINK
withMIRp=apply(interacs,1,function(x) 
	miRtargetsP[miRtargetsP$mature_mirna_id==mirIDs$mature[mirIDs$precursor==x[3]]&miRtargetsP$target_ensembl==x[2],])
temp=interacs[sapply(withMIRp,nrow)>0,]
withMIRp=withMIRp[sapply(withMIRp,nrow)>0]
withMIRp=do.call(rbind,lapply(1:nrow(temp),function(x) 
	cbind(temp[x,c(1,5,3)],withMIRp[[x]][,c(3,1,7,9)])))
colnames(withMIRp)[2]="pam50Symbol"
i=paste(withMIRp$pam50Symbol,withMIRp$predictor)
temp=lapply(unique(i),function(x) withMIRp[i==x,])
withMIRp=do.call(rbind,lapply(temp,function(x)
	cbind(paste(unique(x$subtype),collapse=', '),x[1,2:4],
	      paste(unique(x$experiment),collapse=', '),
	      paste(unique(x$database),collapse=', '),
	      paste(unique(x$pubmed_id),collapse=', '))))
#repeat for validated targets
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
write.table(TFtargets,"TFtargets.tsv",sep='\t',quote=F,row.names=F)

withTF=do.call(rbind,lapply(unique(interacs$pam50Symbol),function(x) 
	interacs[interacs$pam50Symbol==x,][which(interacs$predictorSymbol[interacs$pam50Symbol==x]%in%tfs$TF[tfs$target==x]),]))

######################################################################
#######all 			    
selected0=table(substr(interacs$predictor,1,1))
validated=c(nrow(withCpG),nrow(withTF),nrow(withMIR))
fisher.test(cbind(validated,selected0),alternative="greater")
#p-value < 2.2e-16
count=c(selected,validated)
validated=c(rep("n",3),rep("y",3))
omics=c("CpG","transcript","miRNA","CpG","transcript","miRNA")
selected=as.data.frame(cbind(count,validated,omics),stringsAsFactors=F)
png("ValiOmics.png")
ggplot(selected,aes(fill=omics,y=as.numeric(count),x=validated))+
	geom_bar(position="dodge",stat="identity")+ylab("predictions")+
	scale_y_continuous(trans='log10')+
	scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()


######################################################################
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
revlog_trans <- function(base = exp(1)){
    ## Define the desired transformation.
    trans <- function(x){
                 -log(x, base)
                }
    ## Define the reverse of the desired transformation
    inv <- function(x){
                 base^(-x)
                }
    ## Creates the transformation
    trans_new(paste("revlog-", base, sep = ""),
              trans, ## The transformation function (can be defined using anonymous functions)
              inv,  ## The reverse of the transformation
              log_breaks(base = base), ## default way to define the scale breaks
              domain = c(1e-100, Inf) ## The domain over which the transformation is valued
             )
    }
plotsN=lapply(unique(interacs$predictor),function(x) 
	ggplot(negas[negas$predictor==x,],aes(x=abs(as.numeric(coef))))+
		xlab("coefficient")+ylab("scaled frequency")+
	 	geom_density(aes(group=subtype,color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
	      	scale_x_continuous(trans=revlog_trans(base=10),labels=trans_format("identity",function(x) -x),lim=c(10000,1e-12))+ggtitle(x))
plotsP=lapply(unique(interacs$predictor),function(x)
	ggplot(posi[posi$predictor==x,],aes(x=as.numeric(coef)))+xlab("coefficient")+
	      ylab("scaled frequency")+
	      geom_density(aes(group=subtype,color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
	      scale_x_continuous(trans="log10",lim=c(1e-12,1e6)))
png("OmicsCoeffs.png")
grid.arrange(plotsN[[1]],plotsP[[1]],plotsN[[2]],plotsP[[2]],plotsN[[3]],plotsP[[3]])
dev.off()

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
