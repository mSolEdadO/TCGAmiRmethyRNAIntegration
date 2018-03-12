library(TCGAbiolinks)

#unique patients in the dataset
xC=unique(patients(xprssn))
#as I only have data from the intersection of the three platforms,
#xC serves for all

#files for every patient from each platform
files=lapply(xC,function(x) 
	rbind(getResults(xprssn)[patients(xprssn)==x,c(7,17,27,5,3)],
		getResults(mirna)[patients(mirna)==x,c(8,19,29,6,4)],
		getResults(mthyltn)[patients(mthyltn)==x,c(7,18,27,5,3)]))
#xC[sapply(files, nrow)>6]
#lapply(files[sapply(files, nrow)>6],function(x) table(x[,2:3]))
files=do.call(rbind,files)
files[,4]=paste(files[,4],files[,5],sep='/')
colnames(files)[4]="file"
files=files[,1:4]

#check how many unique transcripts there are per file
#for x in $(cut -f4 archivos.tsv);do cut -f1 $x|sort|uniq|wc -l;done
#added manually to file
files1=read.table("archivos.tsv",sep='\t',header=T)
files1=cbind(substr(files1$cases,1,12),files1)

#plot it
svg("mRNAxPatient.svg")
plot(files1[files1[,4]=="RNA-Seq" &
	!files1[,1]%in%c("TCGA-A7-A13G","TCGA-A7-A13E","TCGA-A7-A0DB") &
	files1[,3]=="Primary solid Tumor",6],
	col="cornflowerblue",type='b',ylim=c(800,1200),
	ylab="found transcripts",xlab="patient index")
lines(files1[files1[,4]=="RNA-Seq" &
	!files1[,1]%in%c("TCGA-A7-A13G","TCGA-A7-A13E","TCGA-A7-A0DB") &
	files1[,3]=="Solid Tissue Normal",6],
	col="coral1",type='b')
legend("bottomright",c("Tumor","Normal tissue"),
	fill=c("cornflowerblue","coral1"),bty='n')
dev.off()
svg("mRNA-TvsN.svg")
plot(files1[files1[,4]=="RNA-Seq" &
	!files1[,1]%in%c("TCGA-A7-A13G","TCGA-A7-A13E","TCGA-A7-A0DB"),
	c(3,6)],col=c("cornflowerblue","coral1"),frame.plot=F,ylim=c(800,1200))
dev.off()
