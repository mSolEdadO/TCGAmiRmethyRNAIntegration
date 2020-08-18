library(ggplot2)
library(gridExtra)
#######################REGULATORS SHARED WITHIN SUTYPES#########################
#regulators shared by BP
temp=lapply(regus,function(x) table(unlist(x)))
#data frame to plot
temp=lapply(temp,function(x) cbind(names(x),x))
temp=data.frame(do.call(rbind,lapply(1:5,function(x) 
 cbind(names(temp)[x],temp[[x]]))),stringsAsFactors=T)
colnames(temp)=c("subtype","name","bp")
temp$regulator=gsub("E","TF",gsub("h","miRNA",gsub("c","CpG",substr(temp$name,1,1))))
temp$regulator=factor(temp$regulator,levels=c("CpG","TF","miRNA"))
temp$bp=as.numeric(as.character(temp$bp))

png("BPpeRegu.png")
 ggplot(temp,aes(y=bp,x=subtype,color=subtype))+geom_boxplot()+
 facet_wrap(~regulator)+xlab("")+ylab("processes by regulator")+
 theme(text=element_text(size=18),axis.text.x=element_blank())
dev.off()

#exclusivity of regulators
#get categories needed
temp$exclu=temp$bp==1
temp$exclu=gsub(TRUE,"exclusive",temp$exclu)
temp$exclu=gsub(FALSE,"shared",temp$exclu)

png("regulatorShared.png")
 ggplot(temp1,aes(y=proportion,fill=regulator,x=subtype))+
 geom_bar(position="fill", stat="identity")+
 facet_wrap(~exclusivity)+xlab("")+ylab("proportion")+
 theme(text=element_text(size=18),axis.text.x=element_text(angle=45))+
 scale_fill_manual(values=gray.colors(5))
dev.off()

#fisher test per omic
temp1=lapply(names(regus),function(x) table(temp[temp$subtype==x,c(4,5)]))
names(temp1)=names(regus)
temp1=lapply(c("CpG","TF","miRNA"),function(x) 
	lapply(temp1,function(y) rbind(y[rownames(y)==x,],colSums(y[rownames(y)!=x,]))))
names(temp1)=c("CpG","TF","miRNA")
#display p.adjusted values in a table subtype vs regulator
matrix(p.adjust(sapply(temp1,function(x) sapply(x,function(y)
 fisher.test(y,alternative="g")$p.val)),"fdr"),ncol=3)#                          CpG          TF      miRNA
matrix(p.adjust(sapply(temp1,function(x) sapply(x,function(y)
 fisher.test(y,alternative="l")$p.val)),"fdr"),ncol=3)

#top shared per omic per subtype
temp=temp[order(temp$bp,decreasing=T),]
temp1=as.matrix(do.call(rbind,lapply(unique(temp$subtype),function(x) 
	do.call(rbind,lapply(unique(temp$regulator),function(y) 
		temp[temp$subtype==x&temp$regulator==y,][1,])))))
#explicit names
temp1[grep("ENSG",temp1[,2]),2]=sapply(temp1[grep("ENSG",temp1[,2]),2],function(x)
	myannot$hgnc_symbol[myannot$ensembl_gene_id==x])
temp1[grep("cg",temp1[,2]),2]=paste("CpG",sapply(temp1[grep("cg",temp1[,2]),2],function(x) 
	unique(methy$V2[methy$V1==x])),sep=":")
temp1=data.frame(temp1)
temp1$bp=as.numeric(as.character(temp1$bp))
#force order CpG->TF->miR
temp1=temp1[order(paste(temp1$subtype,rep(c(2,1,3),5))),]

png("TopRegus-BP.png") 
 ggplot(temp1,aes(x=bp,y=reorder(paste(subtype,name),15:1),fill=subtype))+
 geom_bar(position="dodge", stat="identity")+
 theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
 	 text=element_text(size=18))+ylab("")+xlab("associated processes")+
 geom_text(aes(label=name,hjust="right"),position=position_nudge(x=35))
dev.off()

#jaccard index between BPs
#separate omics
temp=lapply(regus,function(x) lapply(x,function(y) 
	lapply(c("c","E","h"),function(z) y[substr(y,1,1)==z])))
#jaccard per omic
jaccardI=lapply(1:3,function(i) 
	lapply(temp,function(x) 
		sapply(1:length(x),function(y) sapply(1:length(x),function(z)
		 length(intersect(x[[y]][[i]],x[[z]][[i]]))/length(union(x[[y]][[i]],x[[z]][[i]]))))))
#data frame to plot
names(jaccardI)=c("CpG","TF","miRNA")
jaccardI=lapply(jaccardI,function(x) do.call(rbind,lapply(1:5,function(y) 
	cbind(names(x)[y],x[[y]][upper.tri(x[[y]])]))))
jaccardI=data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(names(jaccardI)[x],jaccardI[[x]]))))
colnames(jaccardI)=c("regulator","subtype","jaccardIndex")
jaccardI$jaccardIndex=as.numeric(as.character(jaccardI$jaccardIndex))
jaccardI$regulator=factor(jaccardI$regulator,levels=c("CpG","TF","miRNA"))

png("perSubtypeSharingP.png")
 ggplot(jaccardI,aes(x=jaccardIndex))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~regulator,ncol=1)+theme(text=element_text(size=18))+xlab("jaccard index")
dev.off()

#ks comparison between subtypes per regulator
lapply(unique(jaccardI$regulator),function(z) 
	matrix(round(p.adjust(sapply(names(g),function(x) sapply(names(g),function(y) 
		ks.test(jaccardI$jaccardIndex[jaccardI$regulator==z&jaccardI$subtype==x],
			jaccardI$jaccardIndex[jaccardI$regulator==z&jaccardI$subtype==y])$p.val)),"fdr"),4),ncol=5))
#ks comparison between regulators per subtype
jaccardI=lapply(names(g),function(x) jaccardI[jaccardI$subtype==x,])
lapply(jaccardI,function(x) 
	matrix(round(p.adjust(sapply(unique(x$regulator),function(y) 
		sapply(unique(x$regulator),function(z) 
		ks.test(x$jaccardIndex[x$regulator==y],
			x$jaccardIndex[x$regulator==z])$p.val)),"fdr"),4),ncol=3))

#######################REGULATORS SHARED BETWEEN SUTYPES#########################
#get BPs enriched in several nets
i=names(which(table(unlist(lapply(regus,names)))>1))
temp=lapply(i,function(x) sapply(regus,function(y) y[names(y)==x]))
#filter out subtypes with no regulator for the BP
temp=lapply(temp,function(y) y[sapply(y,function(x) length(unlist(x)))>0])
#how many BPs are shared between subtypes
temp1=table(sapply(temp,function(x) 
	paste(sapply(strsplit(names(x),".",fixed=T),function(y) y[1]),collapse='&')))
#data frame to plot
temp1=data.frame(cbind(names(temp1),temp1))
temp1$temp1=as.numeric(as.character(temp1$temp1))
temp1=temp1[order(temp1$temp1,decreasing=T),]
temp1$V1=factor(temp1$V1,levels=as.character(temp1$V1))

png("BPonSubt.png")
 ggplot(temp1,aes(x=temp1,y=V1))+geom_bar(position="dodge", stat="identity")+
 ylab("")+xlab("count")+theme(text=element_text(size=18))
dev.off()

#separate per omic
temp=lapply(c("c","E","h"),function(z) 
	lapply(temp,function(y) lapply(y,function(x)
	 unlist(x)[substr(unlist(x),1,1)==z])))
#jaccard index 
jaccardI=lapply(1:3,function(i) lapply(temp[[i]],function(z) 
	lapply(1:(length(z)-1),function(x) lapply((x+1):length(z),function(y) 
		rbind(paste(sapply(strsplit(names(z),".",fixed=T),function(a) a[1])[c(x,y)],
			collapse="-"),
			length(intersect(z[[x]],z[[y]]))/length(union(z[[x]],z[[y]])))))))
#data frame to plot
jaccardI=lapply(1:3,function(z) lapply(1:length(i),function(y) 
	lapply(jaccardI[[z]][[y]],function(x) rbind(i[y],x))))
jaccardI=lapply(jaccardI,function(x) t(matrix(unlist(x),nrow=3)))
jaccardI=data.frame(do.call(rbind,lapply(1:3,function(x) 
	cbind(c("CpG","TF","miRNA")[x],jaccardI[[x]]))))
colnames(jaccardI)=c("regulator","bp","pair","jaccardIndex")
jaccardI$jaccardIndex=as.numeric(as.character(jaccardI$jaccardIndex))
#NA arise from comparing 2 empty sets
jaccardI=jaccardI[!is.na(jaccardI$jaccardIndex),]
jaccardI$regulator=factor(jaccardI$regulator,levels=c("CpG","TF","miRNA"))
#trick to plot
jaccardI=jaccardI[order(jaccardI$jaccardIndex,decreasing=T),]
jaccardI$pos=1
jaccardI$pos[jaccardI$regulator=="CpG"]=1:sum(jaccardI$regulator=="CpG")
jaccardI$pos[jaccardI$regulator=="TF"]=1:sum(jaccardI$regulator=="TF")
jaccardI$pos[jaccardI$regulator=="miRNA"]=1:sum(jaccardI$regulator=="miRNA")

cols=c("red4","red2","coral","magenta","darkolivegreen4","darkkhaki","magenta3",
	"green","magenta4","deepskyblue")
#taken from http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
get_legend<-function(myggplot){
   tmp <- ggplot_gtable(ggplot_build(myggplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)}
p=ggplot(jaccardI,aes(x=pos,y=jaccardIndex,labels=bp))+
 	 geom_point(aes(color=pair))+scale_color_manual(values=cols)+
 	 theme(legend.title = element_blank(),legend.position='bottom')
legend=get_legend(p)
plots=lapply(levels(jaccardI$regulator),function(x)
	ggplot(jaccardI[jaccardI$regulator==x,],aes(x=pos,y=jaccardIndex,labels=bp))+
 	 geom_point(aes(color=pair))+geom_text(aes(
 	 	label=ifelse(jaccardIndex[regulator==x]>quantile(jaccardIndex[regulator==x],
 	 		.99),as.character(bp),'')),check_overlap=T,hjust=-0.01)+
 	 ggtitle(x)+xlab("biological process")+ylab("Jaccard Index")+
 	 theme(text=element_text(size=18),axis.text.x=element_blank(),
 	 	axis.ticks.x=element_blank(),legend.position='n',
 	 	legend.title=element_blank())+
 	 scale_color_manual(values=cols))
png("accrossSubtyBP.png",width=1500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],legend,
	ncol=3,nrow=2,layout_matrix=rbind(c(1,2,3),c(4,4)),
	widths=c(3,3,3),heights=c(5,0.4))
dev.off()

#duplicate values for each subtype in the pair	
temp=sapply(strsplit(as.character(jaccardI$pair),"-"),function(x) x[2])	
jaccardI$subtype=sapply(strsplit(as.character(jaccardI$pair),"-"),
	function(x) x[1])
jaccardI=rbind(jaccardI,jaccardI)
jaccardI$subtype[(nrow(jaccardI)-length(temp)+1):nrow(jaccardI)]=temp

#ks comparison between subtypes per regulator
lapply(unique(jaccardI$regulator),function(z) 
	matrix(round(p.adjust(sapply(names(g),function(x) sapply(names(g),function(y) 
		ks.test(jaccardI$jaccardIndex[jaccardI$regulator==z&jaccardI$subtype==x],
			jaccardI$jaccardIndex[jaccardI$regulator==z&jaccardI$subtype==y])$p.val)),"fdr"),4),ncol=5))
#ks comparison between regulators per subtype
jaccardI=lapply(names(g),function(x) jaccardI[jaccardI$subtype==x,])
lapply(jaccardI,function(x) 
	matrix(round(p.adjust(sapply(unique(x$regulator),function(y) 
		sapply(unique(x$regulator),function(z) 
		ks.test(x$jaccardIndex[x$regulator==y],
			x$jaccardIndex[x$regulator==z])$p.val)),"fdr"),4),ncol=3))

#######################OVER/SUB REPRESENTATION#########################
#BPs with/without CpGs
temp=rbind(sapply(regus,function(x) 
	sum(sapply(x,function(y) sum(substr(y,1,1)=="c"))==0)),
	sapply(regus,function(x) sum(sapply(x,function(y) sum(substr(y,1,1)=="c"))>0)))
#test enrichment
fisher.test(cbind(temp[,5],rowSums(temp[,1:4])),alternative="g")

#subrepresentation of miRNAs in BrCa BPs
temp=sapply(regus,function(x) 
	table(unlist(sapply(x,function(y) substr(y,1,1)=="h"))))
#      Basal  Her2 LumA  LumB normal
#FALSE  9097 13970 3292 14370  10790
#TRUE   6134  6708 2018  7731   7957
fisher.test(cbind(rowSums(temp[,1:4]),temp[,5]),alternative="g")
#there are more non-miRNAs in subtypes