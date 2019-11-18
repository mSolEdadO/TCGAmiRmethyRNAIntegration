library(RColorBrewer)
library(gplots)
library(gridExtra)
library(igraph)
library(biomaRt)

#Eval subtype differences
qualy=read.table("modelQuality.tsv",header=T)
boxes=ggplot(qualy,aes(x=subtype,y=RMSE,color=subtype))+geom_boxplot()+ylab("testing RMSE")+
      geom_signif(test='ks.test',
                  comparisons=list(c("LumA","Basal"),c("LumA","LumB"),c("normal","Basal"),c("normal","Her2"),c("normal","LumB")),
                  map_signif_level=T,y_position=c(70000,75000,90000,170000,80000))
densi=ggplot(qualy,aes(x=RMSE))+
       geom_density(aes(group=subtype,color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
       ylab("scaled frequency")+xlab("testing RMSE")+scale_x_continuous(trans='log10')
png("RMSE.png")
 grid.arrange(boxes,densi,nrow=2)
dev.off()
#ks p-value=probability of a test statistic >= than the one observed if samples came from the same distribution
boxes=ggplot(qualy,aes(x=subtype,y=predictors,color=subtype))+geom_boxplot()+
      geom_signif(test='ks.test',
                  comparisons=list(c("LumA","Basal"),c("LumA","Her2"),c("LumA","LumB"),c("LumA","normal"),c("Her2","Basal")),
                  y_position=c(5500,5100,5600,6000,1500),map_signif_level=T)
densi=ggplot(qualy,aes(x=predictors))+
      geom_density(aes(group=subtype,color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
      ylab("scaled frequency")+scale_x_continuous(trans='log10')
png("predisDistri.png")
  grid.arrange(boxes,densi,nrow=2)
dev.off()

g=graph.data.frame(qualy[,1:2],directed=F)
E(g)$weight=qualy$predictors
M=g[unique(qualy$pam50),unique(qualy$subtype)]
M=M[order(match(rownames(M),pam50$ensembl_gene_id)),]
rownames(M)=pam50$hgnc_symbol
png("PrediNum.png")
heatmap.2(as.matrix(M),col=rev(heat.colors(40)),trace="none",colRow=brewer.pal(n=4,name="Set2")[pam50$class],scale="r",key=F,Colv=F,Rowv=F,dendrogram="none",srtCol=45)
legend("topright",fill=brewer.pal(n=4,name="Set2"),legend=levels(pam50$class),bty="n",border="white")
dev.off()
png("RMSEvsNumPredi.png")
 ggplot(qualy,aes(x=RMSE,y=predictors))+geom_point()+ylab("selected predictors")
dev.off()

#eval omic differences
docus=list.files("data",full.names=T)
data=lapply(docus,read.table,skip=1)
ocus=gsub("data/","",docus)
docus=do.call(rbind,lapply(1:250,function(x) cbind(docus[x],data[[x]])))
docus=cbind(t(sapply(strsplit(as.character(docus[,1]),".",fixed=T),function(x) cbind(x[2],x[1]))),docus)
colnames(docus)=c("subtype","pam50","trash","predictor","coef")
docus$trash=NULL
omics=as.data.frame(table(substr(unique(interacs$predictor),1,1)))
omics[,1]=c("CpG","transcript","miRNA")
input=omics
input[,2]=c(384575,16475,433)  
omics=rbind(omics,input)
omics$selected=c(rep("output",3),rep("input",3))
png("SlctdOmic.png")
 ggplot(omics, aes(fill=Var1, y=as.numeric(Freq), x=selected)) +
        geom_bar(position="fill", stat="identity")+ylab("proportion")+xlab("")+
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()
omics1=table(omics[,c(1,3)])
omics1[,2]=omics$Freq[c(1,3,2)]
omics1[,1]=omics$Freq[c(4,6,5)]
omics1=omics1[,2:1]
chisq.test(omics1)
#X-squared = 19692, df = 2, p-value < 2.2e-16
fisher.test(rbind(omics1[3,],colSums(omics1[1:2,])),alternative="greater")
fisher.test(rbind(omics1[2,],colSums(omics1[c(1,3),])),alternative="greater")
fisher.test(rbind(omics1[1,],colSums(omics1[2:3,])),alternative="less")

coefs1=as.matrix(docus)
coefs1[mirs,3]="miRNA"
coefs1[cpgs,3]="CpG"
coefs1[transcri,3]="transcript"
coefs1=as.data.frame(coefs1)
png("coeffs.png")
 boxes=ggplot(coefs1,aes(y=as.numeric(as.character(coef)),color=predictor))+geom_boxplot()+ scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ylab("regression coefficient")
 densi=ggplot(coefs1,aes(x=as.numeric(as.character(coef))))+geom_density(aes(group=predictor,fill=predictor,color=predictor,y=..scaled..),alpha=0.3)+ scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ylab("scaled frequency")+xlim(-0.5,0.5)+xlab("regression coefficient")
 grid.arrange(boxes,densi,nrow=2)
dev.off()
ks.test(docus$coef[cpgs],docus$coef[transcri])
#D = 0.26268, p-value < 2.2e-16
t.test(docus$coef[cpgs],docus$coef[transcri])
#t = 2.0874, df = 19695, p-value = 0.03687
ks.test(docus$coef[cpgs],docus$coef[mirs])
#D = 0.30373, p-value < 2.2e-16
t.test(docus$coef[cpgs],docus$coef[mirs])
#t = 2.21, df = 19861, p-value = 0.02712
ks.test(docus$coef[transcri],docus$coef[mirs])
#D = 0.29293, p-value < 2.2e-16
t.test(docus$coef[transcri],docus$coef[mirs])
#t = 1.4305, df = 918.89, p-value = 0.1529
png("PrediNum1.png")
 heatmap.2(t(table(coefs1[,c(1,3)])),col=rev(heat.colors(10)),scale="n",Rowv=F,Colv=F,trace="none",dendrogram="none",srtCol=45,margins=c(5,10))
dev.off()


mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",host="http://jan2019.archive.ensembl.org")
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=95)
myannot=getBM(
 attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","external_gene_name","entrezgene"), 
 mart=mart)
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
docus=docus[order(docus$pam50),]
temp=table(docus$pam50)
lala=sapply(1:length(temp),function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(temp)[x]],temp[x]))
docus=cbind(docus,unlist(lala))
docus=docus[order(docus$predictor),]
transcri=grep("^E",docus$predictor,perl=T)
temp=table(as.character(docus$predictor[transcri]))
lala=unlist(sapply(1:length(temp),function(x) rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(temp)[x]],temp[x])))
docus=as.matrix(cbind(docus,docus$predictor))
docus[transcri,6]=lala
docus=as.data.frame(docus[order(docus[,2]),])
colnames(docus)[5:6]=c("PAM50Symbol","predictorSymbol")
write.table(docus,"slctdPrdis.tsv",sep='\t',quote=F,row.names=F)
