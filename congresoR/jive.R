library(r.jive)
library(ggsci)#pa repetir la paleta de colores de MFA.R
load("subset/subset.Rda")

#jive necesita una lista con las distintas omicas
omics=list(t(methy),t(rna),t(mirna))
names(omics)=c("methy","rna","mirna")

JIVEomics=jive(omics)
#Final joint rank: 3 , final individual ranks: 72 45 2 
pdf("JIVEomics.pdf")
showVarExplained(JIVEomics)  
showPCA(JIVEomics,n_joint=2,Colors=pal_lancet()(5)[design$subtype],pch=19)
legend("top",horiz=T,legend=levels(subtipos),fill=pal_lancet()(5),bty='n',border=NA)
showPCA(JIVEomics,n_joint=1,n_indiv=c(1,1,1),Colors=pal_lancet()(5)[design$subtype],pch=19)
dev.off()

#pdf("JIVEheatmaps.pdf",height=465,width=705)
#showHeatmaps(JIVEomics)  
#dev.off()

save(JIVEomics,file="useR/jive.Rda")

