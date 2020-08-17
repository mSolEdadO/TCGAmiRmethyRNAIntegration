library(igraph)
library(ggplot2)
library(data.table)

#########################DEGREE#########################
bpG=lapply(BPenriched,function(x) graph.data.frame(x[,c(1,3)],directed=F))

d=lapply(bpG,degree)
d=lapply(d,function(x) x[order(x,decreasing=T)])
#data frame to plot
d=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(d)[x],names(d[[x]]),d[[x]]))))
colnames(d)=c("subtype","name","degree")
#make omics explicit
d$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(d$name,1,1))))
d$omic[d$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"
d$omic=factor(d$omic,levels=c("CpG","TF","miRNA","transcript"))
d$degree=as.numeric(as.character(d$degree))

png("BP.degree.png")
ggplot(d,aes(x=degree))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~omic)+scale_x_continuous(trans="log10")+
 theme(text=element_text(size=18))
dev.off()

#ks comparison between subtypes per regulator
lapply(unique(d$omic),function(z) 
	matrix(round(p.adjust(sapply(names(g),function(x) sapply(names(g),function(y) 
		ks.test(d$degree[d$omic==z&d$subtype==x],
			d$degree[d$omic==z&d$subtype==y])$p.val)),"fdr"),4),ncol=5))
#ks comparison between regulators per subtype
d=lapply(names(g),function(x) d[d$subtype==x,])
lapply(d,function(x) 
	matrix(round(p.adjust(sapply(unique(x$omic),function(y) 
		sapply(unique(x$omic),function(z) 
		ks.test(x$degree[x$omic==y],x$degree[x$omic==z])$p.val)),"fdr"),4),ncol=4))

#########################BETWEENNESS#########################
bet=lapply(bpG,betweenness)
#all over again
bet=lapply(bet,function(x) x[order(x,decreasing=T)])
bet=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(bet)[x],names(bet[[x]]),bet[[x]]))))
colnames(bet)=c("subtype","name","betweenness")
bet$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(bet$name,1,1))))
bet$omic[bet$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"
bet$omic=factor(bet$omic,levels=c("CpG","TF","miRNA","transcript"))
bet$betweenness=as.numeric(as.character(bet$betweenness))

png("BP.betweenness.png")
ggplot(bet,aes(x=betweenness))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~omic)+scale_x_continuous(trans="log10")+
 theme(text=element_text(size=18))
dev.off()

#extra for alternative plot accounting for all the zeros
bet$logRounded=round(log(bet$betweenness,base=10),digits=0)
bet$categories=paste("1e-",bet$logRounded,sep="")
bet$categories[bet$logRounded<0]="<1"
bet$categories[bet$logRounded==-Inf]="0"

temp=lapply(as.character(unique(bet$subtype)),function(x) 
	sapply(as.character(unique(bet$omic)),function(y) 
		cbind(x,y,table(bet$categories[bet$subtype==x&bet$omic==y]))))
temp=data.frame(do.call(rbind,lapply(temp,function(x) 
	do.call(rbind,lapply(x,function(y) cbind(rownames(y),y))))))
colnames(temp)=c("betweenness","subtype","omic","frequency")
temp$frequency=as.numeric(as.character(temp$frequency))
temp$betweenness=factor(temp$betweenness,
	levels=c(0,"<1",levels(temp$betweenness)[3:11]))
temp$omic=factor(temp$omic,levels=c("CpG","TF","miRNA","transcript"))

png("BP.betweenness.alt.png")
 ggplot(temp,aes(x=betweenness,y=frequency,fill=subtype))+
 geom_bar(position="dodge", stat="identity")+scale_y_continuous(trans="log10")+
 facet_wrap(~omic)+theme(text=element_text(size=18),
 	axis.text.x=element_text(angle=45))
dev.off()

lapply(unique(bet$omic),function(z) 
	matrix(round(p.adjust(sapply(names(g),function(x) sapply(names(g),function(y) 
		ks.test(bet$betweenness[bet$omic==z&bet$subtype==x],
			bet$betweenness[bet$omic==z&bet$subtype==y])$p.val)),"fdr"),4),ncol=5))

bet=lapply(names(g),function(x) bet[bet$subtype==x,])
lapply(bet,function(x) 
	matrix(round(p.adjust(sapply(unique(x$omic),function(y) 
		sapply(unique(x$omic),function(z) 
		ks.test(x$betweenness[x$omic==y],
			x$betweenness[x$omic==z])$p.val)),"fdr"),4),ncol=4))


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