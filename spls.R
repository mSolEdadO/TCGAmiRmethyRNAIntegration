load("subset.Rda")
library(mixOmics)
#(ggsci)#pa repetir la paleta de colores de MFA.R

#pimp data
data=list(mRNA=t(rna),miRNA=t(mirna),methylation=t(methy))
lumA=lapply(data,function(x) t(x[,1:331]))
design = matrix(1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
diag(design) = 0
#1=full model, expected to find: 
#1) highly correlated features,
#2) high graph density,
#3) small number of densely connected modules 
#4)poor discrimination, 

#wrapper.sgcca=block.spls
lumA.sgcca=block.spls(X=lumA,
					  indY=1,
					  ncomp=
					  keepX=,
					  keepY=,
					  design=design,
					  scheme = "horst",
					  mode="canonical",# X and Y play a symetric role
					  scale = TRUE,
					  max.iter=1000)

ncomp=PCA
tune.spls to test.keepX
https://github.com/singha53/diablo/blob/master/analyses/benchmarking/benchmarking_enrichmentConnectivity.Rmd