library(parallel)
library(FactoMineR)
#library(TCGAbiolinks)
library(factoextra)
load("porSubti.RData")

#data to needed format
mirna=do.call(cbind,lapply(concatenadas,function(x) x[[1]]))
rna=do.call(cbind,lapply(concatenadas,function(x) x[[3]]))
methy=do.call(cbind,lapply(concatenadas,function(x) x[[2]]))
omics=list(methy=t(methy),rna=t(rna),mirna=t(mirna))
methyDesign=methyDesign[methyDesign$barcode%in%rownames(omics$methy),]

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(FactoMineR)})
#PCA per omic  
PCAomics=parLapply(cl, omics,function(x) PCA(scale(x),ncp=5,graph=F))
stopCluster(cl)  
save(PCAomics,methyDesign,file="PCAomics.RData")
#do omics separate subtypes?
pdf("PCAomics.pdf")
fviz_pca_ind(PCAomics$methy,geom="point",habillage=methyDesign$subtype,palette=ggsci::pal_lancet("lanonc")(5),title="Metilación de CpGs",pch=19)
fviz_eig(PCAomics$methy,ncp=20,title="Metilación de CpGs")
fviz_pca_ind(PCAomics$mirna,geom="point",habillage=methyDesign$subtype,palette=ggsci::pal_lancet("lanonc")(5),title="Expresión de miRNAs",pch=19)
fviz_eig(PCAomics$mirna,ncp=20,title="Expresión de miRNAs")
fviz_pca_ind(PCAomics$rna,geom="point",habillage=methyDesign$subtype,palette=ggsci::pal_lancet("lanonc")(5),title="Expresión de mRNAs",pch=19)
fviz_eig(PCAomics$rna,ncp=20,title="Expresión de mRNAs")
dev.off()		   

#which miRs separate subtypes?
mirContrib=as.data.frame(PCAomics$mirna$var$contrib)
pdf("mirPCAcontrib.pdf")
par(mar=c(3,7,2,2))
 mirContrib=mirContrib[order(mirContrib$Dim.1,decreasing=T),]
 barplot(mirContrib$Dim.1[1:50],log='x',las=2,horiz=T,names.arg=rownames(mirContrib)[1:50],xlab="Dim2 contribution",cex.names=0.8)
 mirContrib=mirContrib[order(mirContrib$Dim.2,decreasing=T),]
 barplot(mirContrib$Dim.2[1:50],log='x',las=2,horiz=T,names.arg=rownames(mirContrib)[1:50],xlab="Dim2 contribution",cex.names=0.8)
dev.off()

multiom=lapply(concatenadas,function(x) t(scale(do.call(rbind,x))))
cl <- makeCluster(no_cores)
clusterEvalQ(cl,{library(FactoMineR)})
#MFA per subtype
MFAmultio=parLapply(cl,multiom,function(x) MFA(x,group=c(446,384575,13904),name.group=c("mirna","methy","rna"),ncp=5,graph=F,type=rep("s",3)))
stopCluster(cl)  
#how redundant are the 3 omics(theis first PCs)? 
pdf("MFAmultio.pdf")
 fviz_screeplot(MFAmultio[[1]],title="LumA")
  fviz_mfa_axes(MFAmultio[[1]],repel=T)
  fviz_mfa_var(MFAmultio[[1]],"group")
  fviz_mfa_ind(MFAmultio[[1]],geom="point")
 fviz_screeplot(MFAmultio[[2]],title="Basal")
  fviz_mfa_axes(MFAmultio[[2]],repel=T)
  fviz_mfa_var(MFAmultio[[2]],"group")
  fviz_mfa_ind(MFAmultio[[2]],geom="point")
 fviz_screeplot(MFAmultio[[3]],title="LumB")
  fviz_mfa_axes(MFAmultio[[3]],repel=T)
  fviz_mfa_var(MFAmultio[[3]],"group")
  fviz_mfa_ind(MFAmultio[[3]],geom="point")
 fviz_screeplot(MFAmultio[[4]],title="Her2")
  fviz_mfa_axes(MFAmultio[[4]],repel=T)
  fviz_mfa_var(MFAmultio[[4]],"group")
  fviz_mfa_ind(MFAmultio[[4]],geom="point")
 fviz_screeplot(MFAmultio[[5]],title="normal")
  fviz_mfa_axes(MFAmultio[[5]],repel=T)
  fviz_mfa_var(MFAmultio[[4]],"group")
  fviz_mfa_ind(MFAmultio[[5]],geom="point")
dev.off()
save(MFAmultio,file="MFAmultio.RData")

#MFAsubti=MFA(scale(do.call(rbind,omics)),group=table(design$subtype),name.group=c("Basal","Her2","LumA","LumB","normal"),ncp=5,graph=F)       
#se lo come la falta de señal de mirna
#save(PCAomics,MFAsubti,file="exploratorio.Rda")

#library(r.jive)
#JIVEomics=jive(omics)
#save(PCAomics,MFAsubti,MFAmultio,JIVEomics,file="exploratorio.Rda")

		
