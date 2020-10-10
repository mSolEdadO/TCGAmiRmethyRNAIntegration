library(ggplot2)
library(gridExtra)
library(pbapply)
library(ggrepel)

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
 fisher.test(y,alternative="g")$p.val)),"fdr"),ncol=3)
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

#######################REGULATOR EDGES SHARED BETWEEN SUTYPES#########################
#only keep normal enriched BPs also present in a subtype
temp=lapply(1:5,function(x) GS_GO_BP[[x]][names(GS_GO_BP[[x]])%in%names(regus[[x]])])
names(temp)=names(top)
temp$normal=temp$normal[names(temp$normal)%in%names(regus$normal)]
temp$normal=temp$normal[order(match(names(temp$normal),names(regus$normal)))]
#get node IDs involved in a BP in a subtype net
temp=lapply(1:5,function(y) lapply(1:length(temp[[y]]),function(x) 
	which(V(g[[y]])$name%in%c(myannot$ensembl_gene_id[myannot$entrezgene%in%temp[[y]][[x]]],
		regus[[y]][[x]]))))
#get subgraphs for these nodes
temp=lapply(1:5,function(x) lapply(temp[[x]],function(y) induced_subgraph(g[[x]],y)))
#recover BP names so every element of the lists is identified
names(temp)=names(top)
names(temp$Basal)=names(regus$Basal)
names(temp$Her2)=names(regus$Her2)
names(temp$LumA)=names(regus$LumA)
names(temp$LumB)=names(regus$LumB)
names(temp$normal)=names(regus$normal)
#get BPs enriched in several nets
i=names(which(table(unlist(lapply(temp,names)))>1))

#how many BPs are shared between subtypes
temp1=table(sapply(i,function(x) paste(names(regus)[sapply(regus,function(y) 
	sum(names(y)==x))>0],collapse='&')))
#data frame to plot
temp1=data.frame(cbind(names(temp1),temp1))
temp1$temp1=as.numeric(as.character(temp1$temp1))
temp1=temp1[order(temp1$temp1,decreasing=T),]
temp1$V1=factor(temp1$V1,levels=as.character(temp1$V1))

png("BPonSubt.png")
 ggplot(temp1,aes(x=temp1,y=V1))+geom_bar(position="dodge", stat="identity")+
 ylab("")+xlab("count")+theme(text=element_text(size=18))
dev.off()

#group subgraphs by repeated BP
temp=lapply(i,function(y) lapply(temp,function(x) x[names(x)==y]))
#filter subtypes without the BP
temp=lapply(temp,function(x) x[sapply(x,length)>0])
names(temp)=i
intersectionCounts=function(subty1,subty2){
	u=ecount(graph.union(subty1,subty2))
	#repeated edges
	edges=get.edgelist(graph.intersection(subty1,subty2))
	#repeated nodes
	nodes=as.character(edges)#not unique(nodes)
	#a regulator linked to two transcripts should be counted twice
	#only count regulatory nodes
	nodes=nodes[nodes%in%unlist(regus)]
	#if no regulatory nodes shared
	if(length(nodes)==0){return(cbind(u,"all",0))}
	#counts per omic
	counts=table(substr(nodes,1,1))
return(cbind(u,names(counts),counts))}
#get jaccardIndex input
acrossSubty=pblapply(temp,function(z) 
	do.call(rbind,lapply(1:(length(z)-1),function(x) 
 		do.call(rbind,lapply((x+1):length(z),function(y) 
			cbind(names(z)[x],names(z)[y],
 				intersectionCounts(z[[x]][[1]],z[[y]][[1]])))))))
#slooow
#data frame to plot
acrossSubty=data.frame(do.call(rbind,lapply(1:length(temp1),function(x)
 cbind(names(temp1)[x],temp1[[x]]))))
colnames(acrossSubty)=c("BP","subtype1","subtype2","Eunion","regulator","Eintersection")
acrossSubty$Eunion=as.numeric(as.character(acrossSubty$Eunion))
acrossSubty$Eintersection=as.numeric(as.character(acrossSubty$Eintersection))
acrossSubty$JaccardIndex=acrossSubty$Eintersection/acrossSubty$Eunion
acrossSubty$regulator=gsub("E","TF",gsub("c","CpG",gsub("h","miRNA",
	acrossSubty$regulator)))
acrossSubty$regulator=factor(acrossSubty$regulator,levels=c("CpG","TF","miRNA","all"))
#trick to plot
acrossSubty=acrossSubty[order(acrossSubty$JaccardIndex,decreasing=T),]
acrossSubty$pos=1
acrossSubty$pos[acrossSubty$regulator=="CpG"]=1:sum(acrossSubty$regulator=="CpG")
acrossSubty$pos[acrossSubty$regulator=="TF"]=1:sum(acrossSubty$regulator=="TF")
acrossSubty$pos[acrossSubty$regulator=="miRNA"]=1:sum(acrossSubty$regulator=="miRNA")
acrossSubty$pair=paste(acrossSubty$subtype1,acrossSubty$subtype2,sep='-')

cols=c("red4","red2","coral","magenta","darkolivegreen4","darkkhaki","magenta3",
	"green","magenta4","deepskyblue")
#taken from http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
get_legend<-function(myggplot){
   tmp <- ggplot_gtable(ggplot_build(myggplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)}
p=ggplot(acrossSubty,aes(x=pos,y=JaccardIndex,labels=BP))+
 	 geom_point(aes(color=pair))+scale_color_manual(values=cols)+
 	 theme(legend.title = element_blank(),legend.position='bottom',text=element_text(size=18))+
 	 guides(colour = guide_legend(override.aes = list(size=3)))
legend=get_legend(p)
plots=lapply(c("CpG","TF","miRNA"),function(x) {
	l=max(acrossSubty$pos[acrossSubty$regulator==x]);
	acrossSubty$pos[acrossSubty$regulator=="all"]=(l+1):
				  (l+sum(acrossSubty$regulator=="all"));
	ggplot(acrossSubty[acrossSubty$regulator%in%c(x,"all"),],
		aes(x=pos,y=JaccardIndex,labels=BP))+
 	 geom_point(aes(color=pair))+geom_text_repel(aes(
 	 	label=ifelse(JaccardIndex>JaccardIndex[regulator==x][6],
 	 		as.character(BP),'')),hjust=-0.02,size=5)+ggtitle(x)+
 	 xlab("biological process")+ylab("Jaccard Index")+
 	 theme(text=element_text(size=18),axis.text.x=element_blank(),
 	 	axis.ticks.x=element_blank(),legend.position='n',
 	 	legend.title=element_blank())+
 	 scale_color_manual(values=cols)})

png("accrossSubtyBP.png",width=1600,height=600)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],legend,
	ncol=3,nrow=2,layout_matrix=rbind(c(1,2,3),c(4,4)),
	widths=c(3,3,3),heights=c(5,0.4))
dev.off()

#one row per subtype
temp=acrossSubty
temp$subtype2=NULL
temp1=acrossSubty
temp1$subtype1=NULL
temp=data.frame(rbind(as.matrix(temp),as.matrix(temp1)))
#one row per regulator
temp1=temp[rep(which(temp$regulator=="all"),3),]
temp1$regulator=as.character(sapply(c("CpG","TF","miRNA"),function(x) rep(x,sum(temp$regulator=="all"))))
temp=rbind(temp[temp$regulator!="all",],temp1)

#ks comparison between subtypes per regulator
lapply(c("CpG","TF","miRNA"),function(z) 
	matrix(round(p.adjust(sapply(names(top),function(x) 
		sapply(names(top),function(y) 
		ks.test(temp$JaccardIndex[temp$regulator==z&temp$subtype1==x],
			temp$JaccardIndex[temp$regulator==z&temp$subtype1==y])$p.val)),"fdr"),4),ncol=5))
#ks comparison between regulators per subtype
temp=lapply(names(top),function(x) temp[temp$subtype1==x,])
lapply(temp,function(x) 
	matrix(round(p.adjust(sapply(unique(x$regulator),function(y) 
		sapply(unique(x$regulator),function(z) 
		ks.test(x$JaccardIndex[x$regulator==y],
			x$JaccardIndex[x$regulator==z])$p.val)),"fdr"),4),ncol=3))

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