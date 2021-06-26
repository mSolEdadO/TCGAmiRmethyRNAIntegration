library(data.table)
library(mixOmics)

raw=fread("Her2.mtrx")
raw=as.matrix(raw[,2:ncol(raw)],rownames=raw$V1)
raw=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,function(x) raw[x[1]:x[2],])
normi=fread("Her2.normalized")
normi=as.matrix(normi[,2:ncol(normi)],rownames=normi$V1)
normi=apply(cbind(c(1,393133,410210),c(393132,410209,410813)),1,function(x) t(normi[x[1]:x[2],]))
names(raw)=c("CpGs","transcripts","miRNAs")
names(normi)=c("CpGs","transcripts","miRNAs")

#op1:cca=wrapper.sgcca(raw,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op2:cca.n=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=F)
#op3:cca.n1=wrapper.sgcca(normalized,penalty=c(1,.5,.3),ncomp=5,scale=T)
#op1 & op3 give the same results BUT op1 ends (cca$crit) faster
#is raw data better???what if such covergence depends on penalties???
#selectVar() is the same than loadings != 0, NOT loadings.star != 0

#note: the penalty parameters will need to be tuned
#sparsity parameters are chosen for each of the 10 MCCV iterations
#using an internal 5-fold CV loop: the parameters that minimize the
#prediction error [Tenenhaus2014]

params_searcher=function(subtype,lCpG,ltrnscri,lmir){
	size=nrow(subtype$miRNAs)
	i=sample(1:size,round(size/2))
	data=lapply(subtype,function(y) y[i,])
	results=wrapper.sgcca(data,penalty=c(lCpG,ltrnscri,lmir),scale=T)
	#ncomp=1 for training
	evar=as.data.frame(do.call(rbind,results$explained_variance))
	features=as.data.frame(do.call(rbind,lapply(results$loadings,
		function(x) apply(x,2,function(y) sum(y!=0)))))
	out=cbind(evar,features,results$penalty)
	out$omic=rownames(out)
	colnames(out)[1:3]=c("explained_variance","nfeatures","sparsity")
return(out)}

#for every omic
temp=lapply(pars,function(y)
lapply(1:20,function(x) params_searcher(traindata,y,1,1)))
temp=do.call(rbind,lapply(temp,function(x) do.call(rbind,x)))
colnames(temp)[4]="sparsity"
temp$sparsity=as.character(temp$sparsity)
ggplot(temp[temp$omic=="CpGs",],aes(x=sparsity,
	y=explained_variance,group=sparsity))+geom_boxplot()
ggplot(temp[temp$omic=="CpGs",],aes(x=sparsity,
	y=nfeatures,group=sparsity))+geom_boxplot()


####################what happens with 