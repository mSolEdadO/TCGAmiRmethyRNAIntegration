library(TCGAbiolinks)
library(SummarizedExperiment)
library(doParallel)
library(data.table)

subtype=read.table("subtype.tsv",header=T,sep='\t')
#get the data
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450",
	barcode=subtype$barcode)
GDCdownload(mthyltn)
#GDCdownload will download 844 files. A total of 119.236093 GB
####mthyltn=GDCprepare(mthyltn) dies
#build matrix
files=list.files()
files=sapply(files,function(x) paste(x,"/",list.files(x),sep=""))
methy=do.call(cbind,pbapply::pbsapply(files,function(x) 
	fread(x,select=2,stringsAsFactors=F)))#2 has beta values
#~44m 03s 
colnames(methy)=sapply(strsplit(files,".",fixed=T),
	function(x) x[6])#barcodes
rownames(methy)=fread(files[1],select=1,skip=1)$V1
dim(methy)
#[1] 485577    845
write.table(methy,"methy.tsv",sep='\t',quote=F)

#######filter NA probes######################################
nas=rowSums(is.na(methy))
#keep only those measured at least in half the samples
methy=methy[nas<ncol(methy)*.5,]
dim(methy)
#[1] 395775    844

############filter noise-prone probes##############################
#polymorphisms may affect DNAm measurements
methy=methy[grep("rs",rownames(methy),invert=T),]
dim(methy)
#[1] 395710    844

#get probes annotation
annot=fread(files[1])
annot=annot[annot$probe%in%rownames(methy),]
sum(annot$probe==rownames(methy)) #order is the same
#[1] 395710

#there're 9 males that may drive fake DM due to chr imbalance 
methy=methy[!annot$Chromosome%in%c("chrX","chrY"),]
dim(methy)
#[1] 385940    844
annot=annot[annot$probe%in%rownames(methy),]
#drop probes with ambigous mapping
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
methy=methy[annot$Chromosome!="*",]
dim(methy)
#[1] 383408    844

#######impute missing data########################################
load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy450.RData")
load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methyN_ori.RData")
#load("GDCdata/TCGA-BRCA/harmonized/DNA_Methylation/methy27.RData")
methySubti=lapply(c("Basal","LumA","LumB","Her2"),function(x) 
	methy450[,colnames(methy450)%in%methyDesign[methyDesign[,4]==x,1]])
methySubti$normal=mthyltN
#methySubti$Her2=methy27
sapply(methySubti,dim)
#      Basal   LumA   LumB normal   Her2
#[1,] 396065 396065 396065 395814 396065
#[2,]    137    337    178     96     75
nas=lapply(methySubti,function(x) rowSums(is.na(x)))
methySubti=lapply(1:5,function(x) methySubti[[x]][nas[[x]]<ncol(methySubti[[x]])*0.75,])
sapply(methySubti,dim)
#      Basal   LumA   LumB normal   Her2
#[1,] 395820 395818 395817 395814 395815
#[2,]    137    337    178     96    75
#impute missing values
no_cores=detectCores()-1
registerDoParallel(cores=no_cores)
cl=makeCluster(no_cores,type="FORK")
methySubti=parLapply(cl,methySubti,function(x) impute.knn(x,k=15,maxp=nrow(x))$data)
names(methySubti)=c("Basal","LumA","LumB","Her2","normal")
#Buuren and Groothuis suggested that 15 variables are generally sufficient for imputation
#maxp forces to impute always with knn, even when there's too much NA in columns/rows
stopCluster(cl)
save(methySubti,methyDesign,file="imptMethy.RData")
#######beta to M values########################################

#Beta values have severe heteroscedasticity for highly methylated or unmethylated CpG sites.
#M-value method provides much better performance in terms of detection rate and true positive rate
#for both highly methylated and unmethylated CpG sites
methySubti=lapply(unique(methyDesign$subtype),function(x) methy[,methyDesign$subtype==x])
methySubti=lapply(methySubti,function(x) x[,!duplicated(substr(colnames(x),1,12))])
names(methySubti)=unique(methyDesign$subtype)
Msubti=lapply(methySubti,beta2m)
sapply(Msubti,dim)
#      Basal   LumA   LumB   Her2 normal
#[1,] 384575 384575 384575 384575 384575
#[2,]    135    331    177     75     96

pdf("methyDistr.pdf")
  plot(density(methySubti$Basal),col="blue",frame.plot=F)
	lines(density(methySubti$LumA),col="purple")
	lines(density(methySubti$LumB),col="green")
	lines(density(methySubti$Her2),col="orange")
	lines(density(methySubti$normal),col="red")
  plot(density(Msubti$Basal),col="blue",frame.plot=F)
  	lines(density(Msubti$LumA),col="purple")
	lines(density(Msubti$LumB),col="green")
	lines(density(Msubti$Her2),col="orange")
	lines(density(Msubti$normal),col="red")
  legend("topright",legend=c("Basal","Her2","LumA","LumB","normal"),fill=c("blue","orange","purple","green","red"),bty="n")
dev.off()
#######correct batches#########################################
#subtype to duplicates
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtype[c(which(!subtype$barcode%in%j),
  as.numeric(sapply(which(subtype$barcode%in%j),rep,2))),1:4]
#needed coz names are not equal to expression data but barcodes do
designMethy$samples=unlist(sapply(designMethy$barcode,function(x) colnames(methy)[i==x][1]))
designMethy$samples[designMethy$barcode%in%j]=rev(colnames(methy)[
	which(i%in%j)])
designMethy=designMethy[,order(match(designMethy$samples,
									colnames(methy)))]




library(NOISeq)
noiseqData=readData(data =do.call(cbind,Msubti), factors =methyDesign)
myPCAshared=dat(noiseqData, type = "PCA", norm = T,logtransf=T)
pdf("mPCA.pdf")
 explo.plot(myPCAshared, samples = c(1,2), plottype = "scores", factor = "subtype")
 explo.plot(myPCAshared, samples = c(3,5), plottype = "scores", factor = "plate")
dev.off()

save(methySubti,methyDesign,Msubti,annot,file="ini/imptMethy.RData")
