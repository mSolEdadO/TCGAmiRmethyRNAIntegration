library(parallel)
library(FactoMineR)
#library(TCGAbiolinks)
library(factoextra)
load("subset/subset.Rda")

#subti=sapply(1:5,function(x) rbind(methy[[x]],rna[[x]],mirna[[x]]))
#names(subti)=names(methy)

#add clinical info to matrixes
#dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical") 
#desi=do.call(rbind,sapply(1:5,function(x) cbind(names(subti)[x],colnames(subti[[x]]))))
#dataClin=dataClin[dataClin$bcr_patient_barcode%in%substr(design[,1],1,12),c(5,6,33,41)]
#temp=t(sapply(substr(design[,2],1,12),function(x) dataClin[dataClin$bcr_patient_barcode==x,]))
#design=cbind(desi,temp[,1:3])

omics=list(t(methy),t(rna),t(mirna))
names(omics)=c("methy","rna","mirna")
cl <- makeCluster(7)
clusterEvalQ(cl,{library(FactoMineR)})
#MFA per subtype using 3 omics 
PCAomics=parLapply(cl, omics,function(x) PCA(x,ncp=5,graph=F))
save(PCAomics,"exploratorio.Rda")
pdf("PCAomics.pdf")
fviz_pca_var(PCAomics$methy,geom="point",col.var=subtipos,palette=ggsci::pal_lancet("lanonc")(5),pch=19,title="Methy")
fviz_eig(PCAomics$methy, addlabels = TRUE,title="Methy")
fviz_pca_var(PCAomics$mirna,geom="point",col.var=subtipos,palette=ggsci::pal_lancet("lanonc")(5),pch=19,title="miRNA")
fviz_eig(PCAomics$mirna, addlabels = TRUE,title="miRNA")
fviz_pca_var(PCAomics$rna,geom="point",col.var=subtipos,palette=ggsci::pal_lancet("lanonc")(5),pch=19,title="mRNA")
fviz_eig(PCAomics$rna, addlabels = TRUE,title="mRNA")
dev.off()		   
MFAsubti=MFA(scale(do.call(rbind,omics)),group=table(design$subtype),name.group=c("Basal","Her2","LumA","LumB","normal"),ncp=5,graph=F)       
#se lo come la falta de señal de mirna
save(PCAomics,MFAsubti,file="exploratorio.Rda")

multio=cbind(as.data.frame(t(scale(do.call(rbind,omics)))),as.factor(unlist(design$subtype)))
MFAmultio=MFA(multio,group=c(39201,1669,210,1),name.group=c("methy","rna","mirna","subtypes"),ncp=5,graph=F,num.group.sup=4,type=c(rep("s",3),"n"))
save(PCAomics,MFAsubti,MFAmultio,file="exploratorio.Rda")
pdf("MFAmultio.pdf")
fviz_mfa_axes(MFAmultio,repel=T)

library(r.jive)
JIVEomics=jive(omics)
save(PCAomics,MFAsubti,MFAmultio,JIVEomics,file="exploratorio.Rda")

		   
		   
methy=t(do.call(rbind,lapply(subti,function(x) x[,1:1588])))
rna=t(do.call(rbind,lapply(subti,function(x) x[,1589:397394])))
mirna=t(do.call(rbind,lapply(subti,function(x) x[,397395:411298])))
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

save(MFAsubti,MFAtotal,MFAomics,MFAmultiomics,subti,file="MFA.Rda")
       
#la omicas coinciden? 
pdf("FactoMineR.pdf")      
fviz_mfa_var(MFAsubti$Her2,palette="hue",legend="bottom",title="Her2","group")
fviz_screeplot(MFAsubti$Her2,title="Her2")
fviz_mfa_var(MFAsubti$normal,palette="hue",legend="bottom",title="normal","group")
fviz_mfa_axes(MFAsubti$Her2,legend="bottom",title="Her2",repel=T)
fviz_contrib(MFAsubti$Her2, choice = "quanti.var", axes = 1, top = 30,palette=scales::hue_pal()(4)[2:4])
fviz_mfa_ind(MFAsubti$Her2,geom="point",partial=rownames(subti$Her2)[sample(1:75,3)],title="Individuals Her2")
dev.off()

data=do.call(cbind,sgccda.res$X)
MFA(data,group=c(100,37,100),name.group=c("mRNA","miRNA","methylation"),ncp=5,graph=F)

pdf("FactoMineR237.pdf")
fviz_mfa_var(MFAomics,palette="lancet",legend="bottom","group",repel=T)
fviz_mfa_axes(MFAomics,repel=T,legend="bottom")
fviz_mfa_ind(MFAomics,habillage=temp1$"sgccda.res$Y",geom="point",palette=ggsci::pal_lancet("lanonc")(5))
		   
fviz_contrib(MFAomics,choice="quanti.var",top=30,axes=1,palette=scales::hue_pal()(4)[4:1])
dev.off()
