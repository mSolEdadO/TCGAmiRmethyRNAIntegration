library(data.table)
subtipos=read.table("subtipos.tsv")
transcri=read.table("normiARSyN.tsv")
mir=read.table("normiARSyNmiR.tsv")
methy=fread("Mvals.tsv")

colnames(mir)=gsub(".","-",colnames(mir),fixed=T)
mirSubti=colnames(mir)
mirSubti=cbind(mirSubti,sapply(strsplit(mirSubti,"-"),function(x) x[length(x)]))
mirSubti=cbind(mirSubti,substr(mirSubti[,1],1,12))
colnames(mirSubti)=colnames(subtipos)
colnames(methy)=gsub(".","-",colnames(methy),fixed=T)
methySubti=colnames(methy)
methySubti=cbind(methySubti,sapply(strsplit(methySubti,"-"),function(x) x[length(x)]))
methySubti=cbind(methySubti,substr(methySubti[,1],1,12))
colnames(methySubti)=colnames(methySubti)

intersec=lapply(as.character(unique(subtipos$definition)),function(x) 
	intersect(intersect(as.character(subtipos$patient)[subtipos$definition==x],as.character(mirSubti$patient)[mirSubti$definition==x]),as.character(methySubti$patient)[methySubti$definition==x]))
names(intersec)=as.character(unique(subtipos$definition))
concatenadas=lapply(1:5,function(x) list(
									methy=as.matrix(methy)[,colnames(methy)%in%methySubti$barcode[methySubti$definition==names(intersec)[x]&methySubti$patient%in%intersec[[x]]]],
									transcri=as.matrix(transcri)[,colnames(transcri)%in%subtipos$barcode[subtipos$definition==names(intersec)[x]&transcri$patient%in%intersec[[x]]]],
									mir=as.matrix(mir)[,colnames(mir)%in%mirSubti$barcode[mirSubti$definition==names(intersec)[x]&mirSubti$patient%in%intersec[[x]]]]))
names(concatenadas)=names(intersec)
concatenadas=lapply(concatenadas,function(x) sapply(x,function(y) y[,order(substr(colnames(y),1,12))]))
sapply(concatenadas,function(x) ncol(x[[1]]))
#  LumA  Basal   LumB   Her2 normal own subtypes
#   331    135    177     75     75 
#  LumA  Basal   LumB   Her2 normal TCGA subtypes
#   395    125    128     45     75 

save(concatenadas,file="porSubti.RData")
#########################################
#lapply(concatenadas,function(x) sapply(x,function(y) summary(as.numeric(y))))
#$normal
#           methy transcri        mir
#Min.    -8.42100 -39100.0  -7370.000
#Max.     8.40100 370900.0 911600.000
#$LumA
#           methy   transcri        mir
#Min.    -8.37300 -1.046e+06 -1.646e+05
#Max.     8.45200  3.383e+06  9.950e+05
#$LumB
#           methy   transcri        mir
#Min.    -8.40200 -1.130e+06 -19390.000
#Max.     8.45600  2.446e+06 937100.000
#$Her2
#          methy   transcri        mir
#Min.    -8.3730 -719500.00 -1.864e+04
#Max.     8.3530 3143000.00  1.068e+06
#$Basal
#          methy   transcri        mir
#Min.    -8.5370 -1.196e+06 -26820.000
#Max.     8.2840  3.294e+06 680800.000

#since omics have different ranges →scale
concatenadas=lapply(concatenadas,function(x) scale(do.call(rbind,x)))
lapply(concatenadas,function(x) summary(as.numeric(x)))
#$normal
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-13.4100  -0.0515  -0.0498   0.0000  -0.0479 413.9000 
#$LumA
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-97.98000  -0.02165  -0.02113   0.00000  -0.02050 298.80000 
#$LumB
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-109.50000   -0.02229   -0.02174    0.00000   -0.02114  243.50000 
#$Her2
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-60.21000  -0.02019  -0.01969   0.00000  -0.01903 242.70000 
#$Basal
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-109.20000   -0.02097   -0.02045    0.00000   -0.01967  281.10000 

