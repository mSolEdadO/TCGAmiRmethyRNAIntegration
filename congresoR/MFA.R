
library(parallel)
library(FactoMineR)
library(TCGAbiolinks)
load("conca/subsetmasdiff.Rda")

subti=sapply(1:5,function(x) rbind(methy[[x]],rna[[x]],mirna[[x]]))
names(subti)=names(methy)

#add clinical info to matrixes
dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical") 
dataClin=dataClin[dataClin$bcr_patient_barcode%in%substr(desi[,2],1,12),c(5,6,33,41)]
desi=do.call(rbind,sapply(1:5,function(x) cbind(names(rna)[x],colnames(rna[[x]]))))
temp=t(sapply(substr(desi[,2],1,12),function(x) dataClin[dataClin$bcr_patient_barcode==x,]))
desi=cbind(desi,temp[,1:3])
subti$LumA=cbind(t(subti$LumA),as.data.frame(apply(desi[desi[,1]=="LumA",c(3:5)],2,function(x) unlist(x))))
subti$Basal=cbind(t(subti$Basal),as.data.frame(apply(desi[desi[,1]=="Basal",c(3:5)],2,function(x) unlist(x))))
subti$LumB=cbind(t(subti$LumB),as.data.frame(apply(desi[desi[,1]=="LumB",c(3:5)],2,function(x) unlist(x))))
subti$Her2=cbind(t(subti$Her2),as.data.frame(apply(desi[desi[,1]=="Her2",c(3:5)],2,function(x) unlist(x))))
subti$normal=cbind(t(subti$normal),as.data.frame(apply(desi[desi[,1]=="normal",c(3:5)],2,function(x) unlist(x))))

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(FactoMineR)})
#MFA per subtype using 3 omics 
MFAsubti=parLapply(cl, subti,function(x) 
	MFA(x,group=c(100,100,35,3),name.group=c("methy","mRNA","miR","clinical"),ncp=5,graph=F,num.group.sup=4,type=c(rep("s",3),"n")))
#pdf("MAFomics.pdf")
#fviz_mfa_axes(MFAsubti[[1]],geom="arrow",legend="bottom",title=names(MFAsubti)[1])#no sirve con sapply
#fviz_mfa_axes(MFAsubti[[2]],geom="arrow",legend="bottom",title=names(MFAsubti)[2])
#fviz_mfa_axes(MFAsubti[[3]],geom="arrow",legend="bottom",title=names(MFAsubti)[3])
#fviz_mfa_axes(MFAsubti[[4]],geom="arrow",legend="bottom",title=names(MFAsubti)[4])
#fviz_mfa_axes(MFAsubti[[5]],geom="arrow",legend="bottom",title=names(MFAsubti)[5])
#dev.off()

temp=cbind(scale(do.call(rbind,subti)[,1:235]),do.call(rbind,subti)[,236:238])
MFAtotal=MFA(temp,group=c(rep(100,3),3),name.group=c("methy","mRNA","miR","clinical"),ncp=5,graph=F,num.group.sup=4,type=c(rep("s",3),"n"))       
      
rna=do.call(cbind,rna)
mirna=do.call(cbind,mirna)
methy=do.call(cbind,methy)
omics=list(methy,rna,mirna)
names(omics)=c("methy","rna","mirna")
#MFA per omic, considering BRCA subtypes
MFAomics=parLapply(cl, omics,function(x) 
	MFA(x,group=c(331,135,177,75,75),name.group=c("LumA","Basal","LumB","Her2","normal"),ncp=5,graph=F))
pdf("MAFsubti.pdf")
fviz_mfa_var(MFAomics[[1]],geom="point",palette="lancet",legend="bottom",title=names(MFAomics)[1])
fviz_mfa_var(MFAomics[[2]],geom="point",palette="lancet",legend="bottom",title=names(MFAomics)[2])
fviz_mfa_var(MFAomics[[3]],geom="point",palette="lancet",legend="bottom",title=names(MFAomics)[3])
dev.off()
stopCluster(cl)
MFAmultiomics=MFA(scale(do.call(rbind,omics)),group=c(331,135,177,75,75),name.group=c("LumA","Basal","LumB","Her2","normal"),ncp=5,graph=F)
#fviz_mfa_var(MFAmultiomics,geom="point",palette="lancet")#se ve igual que el de miR
       
#la omicas coinciden? 
pdf("FactoMineR.pdf")      
fviz_mfa_var(MFAsubti$Her2,palette="hue",legend="bottom",title="Her2","group")
fviz_screeplot(MFAsubti$Her2,title="Her2")
fviz_mfa_var(MFAsubti$normal,palette="hue",legend="bottom",title="normal","group")
fviz_mfa_axes(MFAsubti$Her2,geom="arrow",legend="bottom",title="Her2")
fviz_contrib(MFAsubti$Her2, choice = "quanti.var", axes = 1, top = 30,palette=scales::hue_pal()(4)[2:4])
fviz_mfa_ind(MFAsubti$Her2,geom="point",partial=rownames(subti$Her2)[sample(1:75,3)],title="Individuals Her2")
dev.off()