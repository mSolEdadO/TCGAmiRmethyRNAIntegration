library(TCGAbiolinks)
library(data.table)

subtype=read.table("subtype.tsv",header=T,sep='\t')
#get the data
mthyltn <-  GDCquery(project = "TCGA-BRCA",
	data.category = "DNA Methylation",
	platform="Illumina Human Methylation 450",
	barcode=subtype$samples)
GDCdownload(mthyltn)
#GDCdownload will download 807 files. A total of X GB
####mthyltn=GDCprepare(mthyltn) dies
#build matrix
files=list.files("GDCdata",full.names=T,recursive=T)
methy=do.call(cbind,pbapply::pbsapply(files,function(x) 
	fread(x,select=2,stringsAsFactors=F)))#2 has beta values
#~45m 24s 
colnames(methy)=sapply(strsplit(files,".",fixed=T),
	function(x) x[6])#barcodes
rownames(methy)=fread(files[1],select=1,skip=1)$V1
dim(methy)
#[1] 485577    850
write.table(methy,"methy.tsv",sep='\t',quote=F)
#5751M
#2486M as tar.gz
#######drop probes with too many na########################################
#before noise-prone filter or the process will be slooow

#subtype to duplicates
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtype[c(which(!subtype$samples%in%j),
  as.numeric(sapply(which(subtype$samples%in%j),rep,2))),1:4]
#needed coz names are not equal to expression data but barcodes do
designMethy$barcode=unlist(sapply(designMethy$samples,function(x)
 colnames(methy)[i==x][1]))
designMethy$barcode[designMethy$samples%in%j]=rev(colnames(methy)[
	which(i%in%j)])
designMethy=designMethy[order(match(designMethy$barcode,
	colnames(methy))),]

total=table(subtype$subtype)
#NA per subtype
nas=lapply(names(total),function(x) 
	rowSums(is.na(methy[,designMethy$subtype==x])))
#keep probes with NA in less than 25% of samples of all subtypes
i=unique(unlist(lapply(1:5,function(x) which(nas[[x]]<total[x]*.25))))
methy=methy[i,]
dim(methy)
#[1] 395728    844

############filter noise-prone probes##############################
#get probes annotation
annot=fread(files[1])
colnames(annot)[1]="probe"
annot=annot[annot$probe%in%rownames(methy),]
annot=annot[order(match(annot$probe,rownames(methy))),]
#sum(annot$probe==rownames(methy)) #order is the same

#sex chrs should be dropped if mixed gender
#since this is not the case, keep chrX
#https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
#methy=methy[!annot$Chromosome%in%c("chrX","chrY"),]
#dim(methy)

#drop probes with ambigous mapping
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
methy=methy[annot$Chromosome!="*",]
nrow(methy)
#[1] 393197

#polymorphisms may affect DNAm measurements
methy=methy[grep("rs",rownames(methy),invert=T),]
nrow(methy)
#[1] 393132

#######impute missing data########################################
library(doParallel)
library(impute)

#separate per subtype
methy=lapply(names(total),function(x) 
	methy[,designMethy$subtype==x])
nas=lapply(methy,function(x) rowSums(is.na(x)))
sapply(nas,max)
#[1]  47  19 165  61  49
sapply(nas,function(x) sum(x>0))#cell to estimate
#[1] 12944  8333 24770 14630  5229

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
#alternative
#methy=lapply(methy,function(x) impute.knn(x,k=15,maxp=nrow(x))$data)

write.table(do.call(cbind,methy),"methyNormi.tsv",sep='\t',quote=F)
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
temp$beta=as.numeric(as.character(temp$beta))
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
temp$M=as.numeric(as.character(temp$M))
png("MDistr.png")
ggplot(temp,aes(x=M))+
geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
	alpha=0.3)+theme(text=element_text(size=18))
dev.off()

#######average duplicates########################################
m=do.call(cbind,m)
i=substr(colnames(m),1,19)
j=i[duplicated(i)]
prefi=m[,!i%in%j]
duplis=m[,i%in%j]

temp=do.call(cbind,lapply(j,function(x) 
  rowMeans(duplis[,substr(colnames(duplis),1,19)%in%x])))
#identify samples with barcode 
colnames(temp)=j
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
write.table(final,"methyM.tsv",sep='\t',quote=F)
#2511M when tar.gz
###########check final data#################################
library(NOISeq)

#u're dragging batch effects, use them as covariates when DM
noiseqData = readData(data = final, factor=subtype)
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = T)
#logtransf=F fails & norm=F flattens points
png("PCA_methy_global.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
 factor = "subtype")})
dev.off()
png("PCA_methy_global_stage.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
 factor = "tumor_stage")})
dev.off()
png("PCA_methy_global_race.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
 factor = "race")})
dev.off()
