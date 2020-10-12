library(data.table)
library(ggplot2)
library(fgsea)
library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(HTSanalyzeR)
#load sif files obtained with adj2sif
files=list.files()
files=files[grep("sif",files)]
sif=lapply(files,fread)
names(sif)=gsub(".sif","",files)
#########################TOP MI interactions#########################
sif=lapply(sif,function(x) x[order(x$V2,decreasing=T),])
i=lapply(sif,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
j=unique(i$Basal)
#top 10k of each type of interaction
sif=lapply(1:10,function(x) 
	do.call(rbind,lapply(j,function(y) sif[[x]][which(i[[x]]==y)[1:10000],])))
#1:10 coz there're 5 sif for CpGs & 5 for miR
names(sif)=gsub(".sif","",files)
i=lapply(sif,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
j=j[2:4]
#drop miR-transcript & miR-miR interactions with p-val<10-6
noMIR=lapply(grep("miR",names(sif),invert=T),function(x) 
	sif[[x]][which(i[[x]]%in%j),])
#paste top 10k miR-transcript & miR-miR interactions per subtype
top=lapply(1:5,function(x) rbind(noMIR[[x]],sif[[grep("miR",names(sif))[x]]]))
names(top)=names(sif)[grep("miR",names(sif),invert=T)]
top=lapply(top,function(x) x[!is.na(x$V2),])

#########################PLOT FINAL MI#########################
#get data frame MI per type of interaction per subtype
i=lapply(top,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
i=lapply(1:5,function(x) cbind(i[[x]],top[[x]]$V2))
i=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(top)[x],i[[x]]))))
colnames(i)=c("subtype","type","MI")
#get dataframe with frequencies of rounded MI 
i$MI2=round(as.numeric(as.character(i$MI)),digits=2)
i=lapply(unique(i$type),function(x) i[i$type==x,])
names(i)=sapply(i,function(x) x$type[1])
j=lapply(i,function(y) 
	do.call(rbind,lapply(names(top),function(x) 
		cbind(x,table(y$MI2[y$subtype==x])))))
j=data.frame(do.call(rbind,lapply(1:length(i),function(x) 
	cbind(names(i)[x],rownames(j[[x]]),j[[x]]))))
rownames(j)=NULL
colnames(j)=c("type","MI","subtype","frequency")
j$type=gsub("c","CpG",j$type)#make the omic explicit
j$type=gsub("h","miRNA",j$type)
j$type=gsub("E","transcript",j$type)
j$type=gsub(" ","-",j$type)
j$MI=as.numeric(as.character(j$MI))
j$frequency=as.numeric(as.character(j$frequency))

i=unique(j$type)
#group by freq range
plot1=ggplot(j[j$type%in%i[c(1,2,4)],],aes(x=MI,y=frequency,color=subtype))+
 geom_point()+
 facet_wrap(~type)+
 theme(text=element_text(size=18))
plot2=ggplot(j[j$type%in%i[c(3,5)],],aes(x=MI,y=frequency,color=subtype))+
 geom_point()+
 facet_wrap(~type)+
 theme(text=element_text(size=18))
png("MIfinal.png",width=650)
 gridExtra::grid.arrange(plot1,plot2)
dev.off()
#########################BP OVER-REPRESENTATION#########################
#BPs to search per subtype
load("gseaComplete.RData")#results of enrich_GSEA.R

#NES heatmap
library(igraph)
temp=lapply(fgseaRes,function(x) x[,c(1,5)])
temp=do.call(rbind,lapply(1:4,function(x) cbind(names(top)[x],temp[[x]])))
temp[,2]=Term(temp[,2])
g=graph.data.frame(temp[,1:2],direct=F)temp[,2]=Term(temp[,2])
E(g)$weight=temp[,3]
temp=t(as.matrix(g[unique(temp[,1]),unique(temp[,2])]))
pdf("NES.pdf",height=15)
heatmap.2(temp,Colv=F,scale='n',trace='n',dendrogram='r',col=greenred(100),keysize=0.25,
	key.title="",key.xlab="NES",density.info='n',lmat = rbind(c(0,4),c(2,1),c(0,3)),
	lhei=c(0.5,5,0.2),lwid=c(1,6),srtCol=0,adjCol = c(0.5,1),margins=c(2,25))
dev.off()

gseaBP=lapply(fgseaRes,function(x) Term(x$pathway))
gseaBP$normal=unique(unlist(gseaBP))
GS_GO_BP<-GOGeneSets(species="Hs",ontologies=c("BP"))
names(GS_GO_BP)=Term(names(GS_GO_BP))
GS_GO_BP=lapply(gseaBP,function(x) GS_GO_BP[names(GS_GO_BP)%in%x])
#background entrez
mart=useEnsembl("ensembl",
	dataset="hsapiens_gene_ensembl",
	host="http://apr2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id",
	"hgnc_symbol",
	"chromosome_name",
	"mirbase_id",
	"entrezgene"), mart=mart)
universo=as.character(myannot$entrezgene[!is.na(myannot$entrezgene)])
#entrez from nets
sets=lapply(top,function(x) unique(as.character(c(x$V1,x$V3))))
sets=lapply(sets,function(x) x[grep("ENSG",x)])#only transcripts needed
sets=lapply(sets,function(x) myannot$entrezgene[myannot$ensembl_gene_id%in%x])
sets=lapply(sets,function(x) sapply(x,function(y) y[!is.na(y)]))
sapply(sets,length)
# Basal   Her2   LumA   LumB normal 
#  9621  11905   8632  10323   7394 

#over-representation per subtype
enriquecimientoBP=lapply(1:5,function(x)
	as.data.frame(multiHyperGeoTest(collectionOfGeneSets=GS_GO_BP[[x]],
 				universe= universo, 
 				hits=as.character(sets[[x]]),
 				minGeneSetSize = 15, 
 				pAdjustMethod = "BH")))
enriquecimientoBP=lapply(enriquecimientoBP,function(x) x[x$Adjusted.Pvalue<0.05,])
names(enriquecimientoBP)=names(top)
sapply(enriquecimientoBP,nrow)
# Basal   Her2   LumA   LumB normal 
#   109    119     34    123    159 

#BP enriched sub-sifs
GS_GO_BP=lapply(1:5,function(x) 
	GS_GO_BP[[x]][names(GS_GO_BP[[x]])%in%rownames(enriquecimientoBP[[x]])])
BPenriched=lapply(1:5,function(y) 
		   unique(do.call(rbind,lapply(GS_GO_BP[[y]],function(x) 
			top[[y]][top[[y]]$V1%in%myannot$ensembl_gene_id[myannot$entrezgene%in%x]|
			top[[y]]$V3%in%myannot$ensembl_gene_id[myannot$entrezgene%in%x],]))))
names(BPenriched)=names(top)
names(GS_GO_BP)=names(top)