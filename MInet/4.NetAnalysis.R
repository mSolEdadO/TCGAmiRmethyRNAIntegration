library(igraph)
library(ggplot2)
library(data.table)
library(gridExtra)
#########################DEGREE#########################
bpG=lapply(BPenriched,function(x) graph.data.frame(x[,c(1,3)],directed=F))

d=lapply(bpG,degree)
#get data frame
d=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(d)[x],names(d[[x]]),d[[x]]))))
colnames(d)=c("subtype","name","degree")
#make omics explicit
d$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(d$name,1,1))))
d$omic[d$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"

#data to plot
#get counts
temp=lapply(unique(d$omic),function(y) lapply(names(top),function(x) 
	as.data.frame(table(d$degree[d$subtype==x&d$omic==y]))))
#get frequencies
propor=lapply(temp,function(x) lapply(x,function(y) y$Freq/sum(y$Freq)))
#explicit subtype
temp=lapply(temp,function(x) do.call(rbind,lapply(1:5,function(y) 
	cbind(names(top)[y],x[[y]]))))
#explicit omic
temp=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(unique(d$omic)[x],temp[[x]]))))
temp$Freq=as.numeric(as.character(temp$Freq))
temp$Var1=as.numeric(as.character(temp$Var1))
colnames(temp)[1:2]=c("omic","subtype")
temp$propor=unlist(propor)

png("BP.degree.png")
	ggplot(temp,aes(x=Var1,y=propor,col=subtype))+geom_point()+facet_wrap(~omic)+
	xlab("degree")+ylab("frequency")+scale_x_continuous(trans="log10")+
	scale_y_continuous(trans="log10")+theme(text=element_text(size=18))
dev.off()

#comparison between subtypes per regulator
lapply(unique(d$omic),function(z) 
	matrix(round(p.adjust(sapply(names(top),function(x) sapply(names(top),function(y) 
		wilcox.test(d$degree[d$omic==z&d$subtype==x],
			d$degree[d$omic==z&d$subtype==y],paired=F)$p.val)),"fdr"),4),ncol=5))
#comparison between regulators per subtype
d=lapply(names(top),function(x) d[d$subtype==x,])
lapply(d,function(x) 
	matrix(p.adjust(sapply(unique(x$omic),function(y) 
		sapply(unique(x$omic),function(z) 
		wilcox.test(x$degree[x$omic==y],x$degree[x$omic==z],paired=F)$p.val)),"fdr"),ncol=4))

#########################TRANSITIVITY#########################
cc=lapply(bpG,transitivity,"local")

#all over again
cc=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(top)[x],V(g[[x]])$name,cc[[x]]))))
colnames(cc)=c("subtype","name","transitivity")
cc$omic=substr(cc$name,1,1)
cc$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(cc$name,1,1))))
cc$omic[cc$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"
cc$transitivity=as.numeric(as.character(cc$transitivity))

i=paste(d$subtype,d$name)
j=paste(cc$subtype,cc$name)
cc=cc[order(match(j,i)),]
d=cbind(d,cc$transitivity)
colnames(d)[5]="transitivity"

#########################PATH LENGTH#########################
#get all shortest paths
l=lapply(g,distances)
targets=lapply(1:5,function(x) 
	BPenriched[[x]]$V3[BPenriched[[x]]$V1%in%names(regus[[x]])])
#keep only paths to functional transcripts
l=lapply(1:5,function(x) l[[x]][,colnames(l[[x]])%in%targets[[x]]])

#test distribution's simmilarity across subtypes
sapply(c("c","E","h"),function(z) sapply(1:4,function(x)  
	p.adjust(wilcox.test(as.numeric(l[[x]][substr(rownames(l[[x]]),1,1)==z,]),
		as.numeric(l[[5]][substr(rownames(l[[5]]),1,1)==z,]),paired=F)$p.val,"fdr")))
#test distribution's simmilarity across omics
lapply(1:5,function(z) matrix(p.adjust(sapply(c("c","E","h"),function(x) sapply(c("c","E","h"),function(y) 
	wilcox.test(as.numeric(l[[z]][substr(rownames(l[[z]]),1,1)==x,]),
		as.numeric(l[[z]][substr(rownames(l[[z]]),1,1)==y,]),paired=F)$p.val)),"fdr"),ncol=3))
#compare TF vs transcript
temp=lapply(l,function(x) x[grep("ENSG",rownames(x)),])
p.adjust(sapply(1:5,function(x) 
	wilcox.test(as.numeric(temp[[x]][rownames(temp[[x]])%in%unlist(regus[[x]]),]),
		as.numeric(temp[[x]][!rownames(temp[[x]])%in%unlist(regus[[x]]),]))$p.val),"fdr")
#separate regulators from non-regulators
lr=lapply(1:5,function(x) l[[x]][rownames(l[[x]])%in%unlist(regus[[x]]),])
l=lapply(1:5,function(x) l[[x]][!rownames(l[[x]])%in%unlist(regus[[x]]),])

#data.frame to plot p1
lr=lapply(lr,function(x) do.call(rbind,lapply(c("c","E","h"),function(y) 
	cbind(y,data.frame(table(as.numeric(x[substr(rownames(x),1,1)==y,])))))))
lr=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(top)[x],lr[[x]]))))
colnames(lr)[1:2]=c("subtype","omic")
lr$omic=gsub("c","CpG",gsub("h","miRNA",gsub("E","TF",lr$omic)))
#data.frame to plot p2
l=lapply(l,function(x) data.frame(table(as.numeric(x))))
l=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(top)[x],l[[x]]))))
colnames(l)[1]="subtype"
l$omic="transcript"
l=l[,c(1,4,2,3)]
#bind both
l=rbind(lr,l)

#calc frequency
i=paste(l$subtype,l$omic)
l=lapply(unique(i),function(x) l[i==x,])
propor=unlist(sapply(l,function(x) x$Freq/sum(x$Freq)))
l=do.call(rbind,l)
l$propor=propor
l$omic=factor(l$omic,levels=c("CpG","TF","miRNA","transcript"))

png("BP.l.png") 
 ggplot(l[!is.na(l$Var1),],aes(x=Var1,y=propor,col=subtype))+geom_line()+facet_wrap(~omic)+
 xlab("shortest path length")+ylab("frequency")+theme(text=element_text(size=18))
dev.off()

#########################ADD BP NODES#########################
#change entrez per ensembl IDs
temp=lapply(GS_GO_BP,function(x) lapply(x,function(y) 
	myannot$ensembl_gene_id[myannot$entrezgene%in%y]))
#keep only genes present in the corresponding net
temp=lapply(1:5,function(x) lapply(temp[[x]],function(y) y[y%in%V(g[[x]])$name]))
#add a mock MI, so transcripts are near their function in the picture
temp=lapply(temp,function(x) do.call(rbind,lapply(1:length(x),function(y) 
	cbind(names(x)[y],1,x[[y]]))))
BPenriched=lapply(1:5,function(x) rbind(BPenriched[[x]],temp[[x]]))
names(BPenriched)=names(regus)

g=lapply(BPenriched,function(x) graph.data.frame(x[,c(1,3)],directed=F))
#write sif for cytoscape
lapply(1:5,function(x) write.table(BPenriched[[x]],
	file=paste(names(BPenriched)[x],"temp",sep='.'),quote=F,row.names=F,sep='\t'))
#ENSG of TFs so u can identify them in cytoscape graph
writeLines(myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3],"temp")

#########################NODE PARAMETERS#########################
sapply(g,function(x) length(V(x)))
# Basal   Her2   LumA   LumB normal 
#  8021  11424   5204  10047   9853 
sapply(g,function(x) length(E(x)))
# Basal   Her2   LumA   LumB normal 
# 20356  24602  10548  25162  26367 
sapply(g,function(x) components(x)$no)
# Basal   Her2   LumA   LumB normal 
#     1      1      1      1      1 
sapply(g,function(x) diameter(x))
# Basal   Her2   LumA   LumB normal 
#     9      9      8      9      9 
sapply(g,function(x) transitivity(x,type="global"))