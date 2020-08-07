library(igraph)
library(ggplot2)
#########################ADD BP NODES#########################
#change entrez per ensembl IDs
temp=lapply(GS_GO_BP,function(x) lapply(x,function(y) 
	myannot$ensembl_gene_id[myannot$entrezgene%in%y]))
#keep only genes present in the corresponding net
temp=lapply(1:5,function(x) lapply(temp[[x]],function(y) y[y%in%V(g[[x]])$name]))
#add a mock MI, so transcripts are neat their function in the picture
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
#  9911  12773   6601  11991  13023 
sapply(g,function(x) length(E(x)))
# Basal   Her2   LumA   LumB normal 
# 23842  27233  12422  28853  33569 
sapply(g,function(x) components(x)$no)
# Basal   Her2   LumA   LumB normal 
#     1      1      1      1      1 
sapply(g,function(x) diameter(x))
# Basal   Her2   LumA   LumB normal 
#     8      9      8      9      8 
sapply(g,function(x) transitivity(x,type="global"))
#      Basal        Her2        LumA        LumB      normal 
#0.016844812 0.007697497 0.014794154 0.017727299 0.003401515 

#########################DEGREE#########################
d=lapply(g,degree)
d=lapply(d,function(x) x[order(x,decreasing=T)])
#extract regulators degree
temp=lapply(1:5,function(x) d[[x]][names(d[[x]])%in%unique(unlist(regus[[x]]))])
#build data.frame to plot
temp=lapply(temp,function(x) cbind(substr(names(x),1,1),x))
names(temp)=names(d)
temp=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(temp)[x],temp[[x]]))))
temp$V2=gsub("c","CpG",temp$V2)
temp$V2=gsub("h","miRNA",temp$V2)
temp$V2=gsub("E","TF",temp$V2)
#extract transcripts degree
temp1=lapply(1:5,function(x) d[[x]][!names(d[[x]])%in%unique(unlist(regus[[x]]))])
temp1=lapply(temp1,function(x) x[grep("ENSG",names(x))])
#build data.frame to plot
temp1=lapply(temp1,function(x) cbind(substr(names(x),1,1),x))
names(temp1)=names(d)
temp1=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(temp1)[x],temp1[[x]]))))
temp1$V2=gsub("E","transcript",temp1$V2)
#paste all together
temp=rbind(temp,temp1)
colnames(temp)=c("subtype","omic","degree")
temp$omic=factor(temp$omic,levels=c("CpG","TF","miRNA","transcript"))
temp$degree=as.numeric(as.character(temp$degree))

png("BP.degree.png")
ggplot(temp,aes(x=degree))+
 geom_density(aes(fill=subtype,color=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~omic)+scale_x_continuous(trans="log10")+
 theme(text=element_text(size=18))
dev.off()

#ks comparison between subtypes per regulator
lapply(unique(temp$omic),function(z) sapply(names(g),function(x) 
	sapply(names(g),function(y) 
		p.adjust(ks.test(temp$degree[temp$omic==z&temp$subtype==x],
			temp$degree[temp$omic==z&temp$subtype==y])$p.val,"fdr"))))
#ks comparison between regulators per subtype
temp=lapply(unique(temp$subtype),function(x) temp[temp$subtype==x,])
lapply(temp,function(x) sapply(unique(x$omic),function(y) 
	sapply(unique(x$omic),function(z) 
		p.adjust(ks.test(x$degree[x$omic==y],x$degree[x$omic==z])$p.val,"fdr"))))

#top degree nodes
temp=lapply(temp,function(y) y[names(y)%in%unlist(sapply(temp,function(x) 
	names(which(x>50))))])
names(temp)=names(regus)
temp=temp[sapply(temp,length)>0]
temp=data.frame(do.call(rbind,lapply(1:4,function(x) 
	cbind(names(temp)[x],names(temp[[x]]),temp[[x]]))))

#########################BETWEENNESS#########################
bet=lapply(g,betweenness)
