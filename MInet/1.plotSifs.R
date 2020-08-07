library(data.table)
library(ggplot2)
#load sif files obtained with adj2sif
files=list.files()
files=files[grep("sif",files)]
sif=lapply(files,fread)
names(sif)=gsub(".sif","",files)

#########################PLOT MI DISTRIBUTION#########################
#get the data frame to plot
i=lapply(sif,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
#the first letter of the ID (on V1 & V3) diferentiates the omic
i=lapply(1:5,function(x) cbind(i[[x]],sif[[x]]$V2))
i=data.frame(do.call(rbind,lapply(1:5,function(x) cbind(names(sif)[x],i[[x]]))))
colnames(i)=c("subtype","type","MI")
i$type=gsub("h","miRNA",i$type)#make the omic explicit
#i$type=gsub("c","CpG",i$type)#make the omic explicit
i$type=gsub("E","transcript",i$type)
i$type=gsub(" ","-",i$type)
i$MI=as.numeric(as.character(i$MI))

#png("MItotal.png")
png("MItotal.miR.png",width=650)#invert # if no-miR plot is desired
ggplot(i,aes(x=MI))+
 geom_density(aes(color=subtype,fill=subtype,y=..scaled..),alpha=0.3)+
 facet_wrap(~type)+
 theme(text=element_text(size=18))
dev.off()

table(i[,2:1])
#                       subtype
#type                     Basal    Her2    LumA    LumB  normal
#  miRNA-CpG             475702   67813 2152421   55670 3525434
#  transcript-CpG         51757  315003 1668608   27263  969608
#  transcript-transcript  64908   19474  223124   51461   64355
#density plots hide differences in interaction's number 
#########################MI SCATTER PLOTS#########################
#get dataframe with frequencies of rounded MI 
i$MI2=round(i$MI,digits=2)
i=lapply(unique(i$type),function(x) i[i$type==x,])
names(i)=sapply(i,function(x) x$type[1])
j=lapply(i,function(y) 
	do.call(rbind,lapply(names(sif),function(x) 
		cbind(x,table(y$MI2[y$subtype==x])))))
j=data.frame(do.call(rbind,lapply(1:length(i),function(x) 
	cbind(names(i)[x],rownames(j[[x]]),j[[x]]))))
rownames(j)=NULL
colnames(j)=c("type","MI","subtype","frequency")
j$MI=as.numeric(as.character(j$MI))
j$freq=as.numeric(as.character(j$freq))

#plot1=ggplot(j[j$type=="transcript-transcript",],aes(x=MI,y=frequency,color=subtype))+
#	  geom_point()+
# 	  ggtitle("transcript-transcript")+
#	  theme(text=element_text(size=18))
#plot1=ggplot(j[!j$type=="transcript-transcript",],aes(x=MI,y=frequency,color=subtype))+
#	  geom_point()+
#	  facet_wrap(~type)+
#	  theme(text=element_text(size=18))

#png("MItotal.alt.png")
# gridExtra::grid.arrange(plot2,plot1)
#dev.off()

#invert # if no-miR plot is desired
plot1=ggplot(j[j$type=="miRNA-miRNA",],aes(x=MI,y=frequency,color=subtype))+
	  geom_point()+
 	  ggtitle("miRNA-miRNA")+
	  theme(text=element_text(size=18))
plot2=ggplot(j[j$type=="miRNA-transcript",],aes(x=MI,y=frequency,color=subtype))+
	  geom_point()+
	  ggtitle("miRNA-transcript")+
	  theme(text=element_text(size=18))
png("MItotal.miR.alt.png")
 gridExtra::grid.arrange(plot2,plot1)
dev.off()
