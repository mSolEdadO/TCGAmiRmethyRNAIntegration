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
fviz_screeplot(MFAmultio)
fviz_mfa_axes(MFAmultio,repel=T)
fviz_mfa_var(MFAmultio,"group")
#temp=sample(rownames(MFAmultio$ind$coord),3)
#[1] "TCGA-OL-A66K-01A-11R-A29R-07" "TCGA-D8-A27E-01A-11R-A16F-07"
#[3] "TCGA-AO-A0JI-01A-21R-A056-07"
fviz_mfa_ind(MFAmultio,partial=temp,geom="point")
fviz_contrib(MFAmultio,"group")
fviz_contrib(MFAmultio,"quanti.var",top=50,palette=scales::hue_pal()(3)[3])
dev.off()
		  
		   
#library(r.jive)
#JIVEomics=jive(omics)
#save(PCAomics,MFAsubti,MFAmultio,JIVEomics,file="exploratorio.Rda")

		
