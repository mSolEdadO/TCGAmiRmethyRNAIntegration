library(parallel)
library(FactoMineR)
load("conca/subsetmasdiff.Rda")

subti=sapply(1:5,function(x) rbind(methy[[x]],rna[[x]],mirna[[x]]))
names(subti)=names(methy)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(FactoMineR)})
#MFA per subtype using 3 omics 
MFAresus=parLapply(cl, subti,function(x) 
	MFA(t(x),group=rep(100,3),name.group=c("methy","mRNA","miR"),ncp=10,graph=F))
pdf("MAFomics.pdf",height=6,width=22,pointsize=26)
par(mfrow=c(1,5))
sapply(1:4,function(x) plot(MFAresus[[x]],
	choix="var",#axes me pone los componentes de c/omica, var me pone todos los vectores
	lab.var=F,
	title=names(subti)[x],
	axes=c(1,2))) #sobre los componentes 1 y 2 del pca conjunto
dev.off()

rna=do.call(cbind,rna)
mirna=do.call(cbind,mirna)
methy=do.call(cbind,methy)
omics=list(methy,rna,mirna)
names(omics)=c("methy","rna","mirna")
#MFA per omic, considering BRCA subtypes
MFAomics=parLapply(cl, omics,function(x) 
	MFA(x,group=c(331,135,177,75,75),name.group=c("LumA","Basal","LumB","Her2","normal"),ncp=10,graph=F))
pdf("MAFsubti.pdf",height=6,width=12,pointsize=20)
par(mfrow=c(1,3))
sapply(1:3,function(x) plot(MFAomics[[x]],choix="var",title=names(omics)[x],lab.var=F,axes=c(1,2)))
dev.off()

#conde example
cg00428457=condes(t(subti$Her2),num.var=1)$quanti
cg00428457=cbind(cg00428457,p.adjust(cg00428457[,2],"fdr"))
colnames(cg00428457)[3]="fdr.adj"
ENSG00000004838=condes(t(subti$Her2),num.var=101)$quanti
ENSG00000004838=cbind(ENSG00000004838,p.adjust(ENSG00000004838[,2],"fdr"))
colnames(ENSG00000004838)[3]="fdr.adj"
hsalet7c=condes(t(subti$Her2),num.var=201)$quanti
hsalet7c=cbind(hsalet7c,p.adjust(hsalet7c[,2],"fdr"))
colnames(hsalet7c)[3]="fdr.adj"

#dimdesc example
temp=dimdesc(MFAresus$LumA,axes=1:2)
temp[[1]]$quanti=cbind(temp[[1]]$quanti,p.adjust(temp[[1]]$quanti[,2],"fdr"))
temp[[2]]$quanti=cbind(temp[[2]]$quanti,p.adjust(temp[[2]]$quanti[,2],"fdr"))
temp[[1]]$quanti[abs(temp[[1]]$quanti[,1])>0.5&temp[[1]]$quanti[,3]<0.05,]
temp[[2]]$quanti[abs(temp[[2]]$quanti[,1])>0.5&temp[[2]]$quanti[,3]<0.05,]