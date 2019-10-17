library(data.table)
interacs=fread("slctdPrdctrs0.tsv")
readHDIout=function(file){
 hdi=fread(file)
 stability=as.data.frame(cbind(hdi[,1],rowSums(hdi[,2:101]),hdi[,102]))
 colnames(stability)[c(1,3)]=c("predictor","pval")
 modelo=unlist(strsplit(file,".",fixed=T))
 stability[stability$predictor%in%interacs$predictor[interacs$subtype==modelo[1]&interacs$pam50==modelo[2]],]}

files=list.files()
files=files[grep("pval",files)]
pvals=lapply(files,readHDIout)
files=files[sapply(pvals,nrow)>0]
pvals=pvals[sapply(pvals,nrow)>0]
pvals=do.call(rbind,lapply(1:length(files),function(x) cbind(files[x],pvals[[x]])))
pvals$binomPval=sapply(pvals$stability,function(y) binom.test(x=y,n=100,p=1/401482,alternative="greater",conf.level=0.99)$p.val)
pvals$fdr.q=p.adjust(pvals$binomPval,"fdr")
