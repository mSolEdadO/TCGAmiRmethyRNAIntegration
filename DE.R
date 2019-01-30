library(limma)

rna=cbind(sapply(concatenadas,function(x) x[[3]]))
design=as.data.frame(do.call(rbind,sapply(1:5,function(x) cbind(colnames(concatenadas[[x]][[1]]),names(concatenadas)[x]))))
colnames(design)=c("sample","subtype")
DE.design=model.matrix(~0+design$subtype)
colnames(DE.design)=gsub("design.","",colnames(DE.design))
contr.mtrx=makeContrasts(
	basal_normal=subtypeBasal-subtypenormal,
	her2_normal=subtypeHer2-subtypenormal,
	luma_normal=subtypeLumA-subtypenormal,
	lumb_normal=subtypeLumB-subtypenormal,
levels=DE.design)	
fitSubtype = contrasts.fit(fit, contr.mtrx)
fitSubtype = eBayes(fitSubtype)
DEGlimma=decideTests(fitSubtype,method="separate",p.value=1e-100,adjust.method="fdr")
vennDiagram(DEGlimma)

#contr.mtrx=makeContrasts(
#	cancer_normal=(subtypeBasal-subtypenormal)-(subtypeHer2-subtypenormal)-(subtypeLumA-subtypenormal)-(subtypeLumB-subtypenormal),
#	cancer_basal=(subtypeBasal-subtypeHer2)-(subtypeBasal-subtypeLumA)-(subtypeBasal-subtypeLumB),
#	cancer_LumA=(subtypeLumA-subtypeHer2)-(subtypeLumA-subtypeLumB),
#	cancer_LumB=(subtypeLumB-subtypeHer2),
#	levels=DE.design)
#fitSubtype = contrasts.fit(fit, contr.mtrx)
#fitSubtype = eBayes(fitSubtype)
#DEGlimma=decideTests(fitSubtype,method="separate",p.value=1e-100,adjust.method="fdr")#https://github.com/CSB-IG/tcgarnaseqbc/blob/master/DifGenes.R
#vennDiagram(DEGlimma)
