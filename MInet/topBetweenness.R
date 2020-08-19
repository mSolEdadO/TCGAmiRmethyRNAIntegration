library(igraph)
#load sif files obtained with adj2sif
topG=lapply(top,function(x) graph.data.frame(x[,c(1,3)],directed=F))
btnss=lapply(topG,betweenness)

btnss=lapply(btnss,function(x) x[order(x,decreasing=T)])
btnss=data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(names(btnss)[x],names(btnss[[x]]),btnss[[x]]))))
colnames(btnss)=c("subtype","name","betweenness")
btnss$omic=gsub("E","transcript",gsub("c","CpG",gsub("h","miRNA",substr(btnss$name,1,1))))
btnss$omic[btnss$name%in%myannot$ensembl_gene_id[myannot$hgnc_symbol%in%tfs$V3]]="TF"
btnss$omic=factor(btnss$omic,levels=c("CpG","TF","miRNA","transcript"))
btnss$betweenness=as.numeric(as.character(btnss$betweenness))
btnss$logRounded=round(log(btnss$betweenness,base=10),digits=0)
btnss$categories=paste("1e-",btnss$logRounded,sep="")
btnss$categories[btnss$logRounded<0]="<1"
btnss$categories[btnss$logRounded==-Inf]="0"

temp=lapply(as.character(unique(btnss$subtype)),function(x) 
	sapply(as.character(unique(btnss$omic)),function(y) 
		cbind(x,y,table(btnss$categories[btnss$subtype==x&btnss$omic==y]))))
temp=data.frame(do.call(rbind,lapply(temp,function(x) 
	do.call(rbind,lapply(x,function(y) cbind(rownames(y),y))))))
colnames(temp)=c("betweenness","subtype","omic","frequency")
temp$frequency=as.numeric(as.character(temp$frequency))
temp$betweenness=factor(temp$betweenness,
	levels=c(0,"<1",levels(temp$betweenness)[3:11]))
temp$omic=factor(temp$omic,levels=c("CpG","TF","miRNA","transcript"))

png("topBetweenness.png")
 ggplot(temp,aes(x=betweenness,y=frequency,fill=subtype))+
 geom_bar(position="dodge", stat="identity")+scale_y_continuous(trans="log10")+
 facet_wrap(~omic)+theme(text=element_text(size=18),
 	axis.text.x=element_text(angle=45))
dev.off()

lapply(names(top),function(x) 
	sapply(c("CpG","TF","miRNA"),function(y) 
		ks.test(bet$betweenness[bet$subtype==x&bet$omic==y],
			btnss$betweenness[btnss$subtype==x&btnss$omic==y])$p.val))