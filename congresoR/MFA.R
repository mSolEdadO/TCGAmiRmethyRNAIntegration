
library(parallel)
library(FactoMineR)
library(TCGAbiolinks)
load("conca/subsetmasdiff.Rda")

subti=sapply(1:5,function(x) rbind(methy[[x]],rna[[x]],mirna[[x]]))
names(subti)=names(methy)

#add clinical info to matrixes
dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical") 
dataClin=dataClin[dataClin$bcr_patient_barcode%in%substr(desi[,2],1,12),]
desi=do.call(rbind,sapply(1:5,function(x) cbind(names(rna)[x],colnames(rna[[x]]))))
temp=t(sapply(substr(desi[,2],1,12),function(x) dataClin[dataClin$bcr_patient_barcode==x,]))
desi=cbind(desi,temp)
subti$LumA=rbind(subti$LumA,t(desi[desi[,1]=="LumA",c(3:6)]))
subti$Basal=rbind(subti$Basal,t(desi[desi[,1]=="Basal",c(3:6)]))
subti$LumB=rbind(subti$LumB,t(desi[desi[,1]=="LumB",c(3:6)]))
subti$Her2=rbind(subti$Her2,t(desi[desi[,1]=="Her2",c(3:6)]))

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(FactoMineR)})
#MFA per subtype using 3 omics 
MFAsubti=parLapply(cl, subti,function(x) 
	MFA(t(unlist(x)),group=c(rep(100,3),4),name.group=c("methy","mRNA","miR","clinical"),ncp=10,graph=F,num.group.sup=301:304,type=c(rep("s",3),rep("n",4))))
pdf("MAFomics.pdf",height=6,width=22,pointsize=26)
par(mfrow=c(1,5))
sapply(1:4,function(x) plot(MFAsubti[[x]],
	choix="axes",#axes me pone los componentes de c/omica, var me pone todos los vectores
	lab.var=F,
	title=names(subti)[x],
	axes=c(1,2))) #sobre los componentes 1 y 2 del pca conjunto
dev.off()
       
######mejor usa axes, para exponer lo de pCA de pCAs-----------------------------       

       
       
       
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
       
####CHECA SI LOS cg ASOCIADOS ESTÁN CERCA, CORRESPONDEN A UN SHORE...?       
       
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
temp=dimdesc(MFAsubti$LumA,axes=1:2)
temp[[1]]$quanti=cbind(temp[[1]]$quanti,p.adjust(temp[[1]]$quanti[,2],"fdr"))
temp[[2]]$quanti=cbind(temp[[2]]$quanti,p.adjust(temp[[2]]$quanti[,2],"fdr"))
temp[[1]]$quanti[abs(temp[[1]]$quanti[,1])>0.5&temp[[1]]$quanti[,3]<0.05,]
temp[[2]]$quanti[abs(temp[[2]]$quanti[,1])>0.5&temp[[2]]$quanti[,3]<0.05,]
