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
#get frequencies
temp=lapply(unique(d$omic),function(y) lapply(names(top),function(x) 
	as.data.frame(table(d$degree[d$subtype==x&d$omic==y]))))
#explicit subtype
temp=lapply(temp,function(x) do.call(rbind,lapply(1:5,function(y) 
	cbind(names(top)[y],x[[y]]))))
#explicit omic
temp=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(unique(d$omic)[x],temp[[x]]))))
temp$Freq=as.numeric(as.character(temp$Freq))
temp$Var1=as.numeric(as.character(temp$Var1))
colnames(temp)[1:2]=c("omic","subtype")

png("BP.degree.png")
	ggplot(temp,aes(x=Var1,y=Freq,col=subtype))+geom_point()+facet_wrap(~omic)+
	xlab("degree")+ylab("frequency")+scale_x_continuous(trans="log10")+
	scale_y_continuous(trans="log10")+theme(text=element_text(size=18))
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

#########################TRANSITIVITY#########################
cc=lapply(bpG,transitivity,"local")

#all over again
cc=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(top)[x],V(g[[x]])$name,cc[[x]]))))
colnames(cc)=c("subtype","name","transitivity")
cc$omic=substr(cc$name,1,1)
cc$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(cc$name,1,1))))
cc$omic[cc$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"
cc$transitivity=as.numeric(as.character(cc$transitivity))

temp=lapply(unique(cc$omic),function(y) lapply(names(top),function(x) 
	as.data.frame(table(cc$transitivity[cc$subtype==x&cc$omic==y]))))
#explicit subtype
temp=lapply(temp,function(x) do.call(rbind,lapply(1:5,function(y) 
	cbind(names(top)[y],x[[y]]))))
#explicit omic
temp=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(unique(cc$omic)[x],temp[[x]]))))
temp$Freq=as.numeric(as.character(temp$Freq))
colnames(temp)[1:2]=c("omic","subtype")


lapply(unique(bet$omic),function(z) 
	matrix(round(p.adjust(sapply(names(g),function(x) sapply(names(g),function(y) 
		ks.test(bet$betweenness[bet$omic==z&bet$subtype==x],
			bet$betweenness[bet$omic==z&bet$subtype==y])$p.val)),"fdr"),4),ncol=5))

bet=lapply(names(top),function(x) bet[bet$subtype==x,])
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