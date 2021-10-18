temp=read.table("penalty_search.tsv",sep='\t',header=T)

#######################################DIAGNOSTIC PLOTS
library(ggplot2)
library(gridExtra)

temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
temp=temp[order(temp$penalty),]
png("AVE.png",width=800)
 ggplot(temp,aes(y=AVE,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))
dev.off()
png("nfeatures.png",width=800)
 ggplot(temp,aes(y=nfeatures,x=as.character(round(penalty,2)),
 group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
 theme(text=element_text(size=18),
 axis.text.x = element_text(angle = 45))+
 scale_y_continuous(trans="log10")
dev.off()
#plot median AVE vs meadian nfeatures
omics=levels(temp$omic)
omics=lapply(omics,function(x) temp[temp$omic==x,])
omics=lapply(omics,function(x) as.data.frame(apply(x,2,as.numeric)))
names(omics)=levels(temp$omic)
omics=lapply(omics,function(y) sapply(unique(y$penalty),function(x) 
apply(y[y$penalty==x,],2,median)))#better than mean?
omics=lapply(omics,function(x) as.data.frame(t(x)))
#indi plots or CpGs will determine axis
plots=lapply(1:3,function(x) ggplot(omics[[x]],
	aes(x=nfeatures,y=AVE,col=penalty))+geom_point()+
	ggtitle(names(omics)[x])+theme(text=element_text(size=18))+
	scale_x_continuous(trans="log10")+geom_line())
png("sparsity_search.png")
 grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

#######################################PENALTIES SUGGESTED BY PLOTS
grid=unique(temp$penalty)
#slopes=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
#	y$AVE[x+1]-y$AVE[x]))
slopes1=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
	abs(y$AVE[x+1]-y$AVE[x])/abs(y$nfeatures[x+1]-y$nfeatures[x])))

slopes1$miRNAs[abs(slopes1$miRNAs)=="Inf"]=NA
#slopes1$transcripts[abs(slopes1$transcripts)=="Inf"]=NA
#penalty=sapply(slopes,function(x) grid[which.max(x)+1])
#       CpGs transcripts      miRNAs 
#       0.20        0.02        0.20 
#slopes1
#       CpGs transcripts      miRNAs 
#       0.02        0.02        0.05 