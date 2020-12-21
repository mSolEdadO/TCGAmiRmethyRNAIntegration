library(TCGAbiolinks)
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
#2600M when tar.gz
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
library(doParallel)
library(impute)

#subtype to duplicates
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtype[c(which(!subtype$barcode%in%j),
  as.numeric(sapply(which(subtype$barcode%in%j),rep,2))),1:4]
#needed coz names are not equal to expression data but barcodes do
designMethy$samples=unlist(sapply(designMethy$barcode,function(x)
 colnames(methy)[i==x][1]))
designMethy$samples[designMethy$barcode%in%j]=rev(colnames(methy)[
	which(i%in%j)])
designMethy=designMethy[order(match(designMethy$samples,
	colnames(methy))),]
#separate per subtype
methy=lapply(levels(designMethy$subtype),function(x) 
	methy[,designMethy$subtype==x])
sapply(methy,dim)
#      Basal   LumA   LumB normal   Her2
#[1,] 383408 383408 383408 383408 383408
#[2,]    126     46    425    146    101
nas=lapply(methy,function(x) rowSums(is.na(x)))
sapply(nas,max)
#[1]  69  23 203  75  55
sapply(nas,function(x) sum(x>0))#cell to estimate
#[1] 12655  8196 24494 15291  9498

#impute missing values
no_cores=detectCores()-1
registerDoParallel(cores=no_cores)
cl=makeCluster(no_cores,type="FORK")
#takes a looong time
methy=parLapply(cl,methy,function(x) 
	impute.knn(x,k=15,maxp=nrow(x))$data)#maxp forces impute, even if too many NAs
#Buuren and Groothuis suggested that 15 variables are generally 
#sufficient for imputation
stopCluster(cl)
write.table(do.call(cbind,methy),"methy.tsv",sep='\t',quote=F)
#######beta to M values########################################
#Beta values have severe heteroscedasticity for highly methylated or
#unmethylated CpG sites. M-values provide much better performance in
#terms of detection rate and true positive rate for both highly 
#methylated and unmethylated CpG sites
library(ggplot2)

#check beta distributions
temp=lapply(methy,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(temp[[x]],levels(designMethy$subtype)[x]))))
temp$beta=as.numeric(as.charachter(temp$beta))
colnames(temp)=c("beta","subtype")
png("betaDistr.png")
ggplot(temp,aes(x=beta))+
geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
	alpha=0.3)+theme(text=element_text(size=18))
dev.off()

#transform beta to M as in:
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
mval=function(beta){log2(beta/(1-beta))}
m=pbapply::pblapply(methy,function(x) apply(x,c(1,2),mval))
#check M distributions
temp=lapply(m,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:5,function(x) 
	cbind(temp[[x]],levels(designMethy$subtype)[x]))))
colnames(temp)=c("M","subtype")
temp$M=as.numeric(as.charachter(temp$M))
png("MDistr.png")
ggplot(temp,aes(x=M))+
geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
	alpha=0.3)+theme(text=element_text(size=18))
dev.off()

#######average duplicates########################################
m=do.call(cbind,m)
i=designMethy$samples[designMethy$barcode%in%j]
prefi=m[,!colnames(m)%in%i]
duplis=m[,colnames(m)%in%i]

temp=do.call(cbind,lapply(j,function(x) 
  rowMeans(duplis[,substr(colnames(duplis),1,19)%in%x])))
#identify samples with barcode 
colnames(temp)=j
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
write.table(final,"methyM.tsv",sep='\t',quote=F)
#2511M when tar.gz

#u're dragging batch effects, use them as covariates when DM