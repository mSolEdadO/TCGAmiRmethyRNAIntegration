library(r.jive)
library(ggsci)#pa repetir la paleta de colores de MFA.R
load("subset/subset.Rda")

#jive necesita una lista con las distintas omicas
omics=list(t(methy),t(rna),t(mirna))
names(omics)=c("methy","rna","mirna")

#hay que ordenar las samples para que los colores por subtipo correspondan con el de MFA.R
omics$methy=omics$methy[,order(match(subtipos,names(table(subtipos))))]
omics$rna=omics$rna[,order(match(subtipos,names(table(subtipos))))]
omics$mirna=omics$mirna[,order(match(subtipos,names(table(subtipos))))]
subtipos=subtipos[order(match(subtipos,names(table(subtipos))))]


JIVEomics=jive(omics)
#Final joint rank: 2 , final individual ranks: 6 3 4 
pdf("JIVEomics.pdf")
showVarExplained(JIVEomics)  
showPCA(JIVEomics,n_joint=2,Colors=pal_lancet()(5)[as.factor(subtipos)],pch=17)
legend("top",horiz=T,legend=unique(subtipos),fill=pal_lancet()(5),bty='n',border=NA)
showPCA(JIVEomics,n_joint=1,n_indiv=c(1,1,1),Colors=pal_lancet()(5)[as.factor(subtipos)],pch=17)
dev.off()

#pdf("JIVEheatmaps.pdf",height=465,width=705)
#showHeatmaps(JIVEomics)  
#dev.off()

save(JIVEomics,file="useR/jive.Rda")

