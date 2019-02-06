library(r.jive)
load("useR/MFA.Rda")

#jive necesita una lista con las distintas omicas
methy=t(do.call(cbind,lapply(subti,function(x) x[,1:100])))
rna=t(do.call(rbind,lapply(subti,function(x) x[,101:200])))
mirna=t(do.call(rbind,lapply(subti,function(x) x[,201:235])))
omics=list(methy,rna,mirna)
names(omics)=c("methy","rna","mirna")
subtipos=gsub(".TCGA.+","",colnames(omics$methy))

JIVEomics=jive(omics)
#Final joint rank: 2 , final individual ranks: 6 3 4 
pdf("JIVEomics.pdf")
showVarExplained(JIVEomics)  
showPCA(JIVEomics,n_joint=2,Colors=scales::hue_pal()(5)[as.factor(subtipos)],pch=17)
legend("top",horiz=T,legend=unique(subtipos),fill=scales::hue_pal()(5),bty='n',border=NA)
showPCA(JIVEomics,n_joint=1,n_indiv=c(1,1,1),Colors=scales::hue_pal()(5)[as.factor(subtipos)],pch=17)
dev.off()

#pdf("JIVEheatmaps.pdf",height=465,width=705)
#showHeatmaps(JIVEomics)  
#dev.off()

save(JIVEomics,file="useR/jive.Rda")

