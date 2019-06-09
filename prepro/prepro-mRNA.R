#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
## &&
## Data preparation: 
##      -Quality Control & bias removal
## By Cristobal Fresno - cristobalfresno@gmail.com

library(SummarizedExperiment)
library(TCGAbiolinks)
#library(Venn.Diagram)
library(biomaRt)  
library(NOISeq)
library(edgeR)
library(cqn)
#library(DESeq2)
library(pbcmc)
library("BiocParallel")
library(sva)

#get the data
xprssn <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary solid Tumor",
  workflow.type = "HTSeq - Counts")
#GDCdownload(xprssn)
expr=GDCprepare(xprssn)
xprssN <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = c("Solid Tissue Normal"),
  workflow.type = "HTSeq - Counts")
#GDCdownload(xprssN)
normal=GDCprepare(xprssN)
duplicados=colData(expr)$patient[duplicated(colData(expr)$patient)]

#table to check for batch effects
designExp=colData(expr)
designExp=rbind(designExp,colData(normal))
#expression matrix: transcripts per sample with the transcript count in each cell
expr=assay(expr)
expr=cbind(expr,assay(normal))


#check the intersection between normal & tumor samples
venn.diagram(x = list(A=designExp$patient[designExp$definition=="Primary solid Tumor"],
  B=designExp$patient[designExp$definition=="Solid Tissue Normal"]), filename = "Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))

#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", "percentage_gene_gc_content", "gene_biotype",
  "start_position","end_position","hgnc_id","hgnc_symbol"),filters = "ensembl_gene_id", 
  values=rownames(expr), mart=mart)
myannot$lenght=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=expr[rownames(expr)%in%myannot$ensembl_gene_id,]
#19172 transcripts

#sum duplicated probes = 2 probes mapping to the same hgnc_id
myannot$hgnc_id[duplicated(myannot$hgnc_id)]
#[1] "HGNC:30046"
myannot[myannot$hgnc_id=="HGNC:30046",]
#      ensembl_gene_id percentage_gene_gc_content   gene_biotype start_position
#41743 ENSG00000254093                      43.17 protein_coding       10764963
#44594 ENSG00000258724                      43.86 protein_coding       10725399
#      end_position    hgnc_id hgnc_symbol lenght
#41743     10839884 HGNC:30046       PINX1  74921
#44594     10839847 HGNC:30046       PINX1 114448
exprots_hgnc[rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"],][1,]=colSums(exprots_hgnc[rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"],])
exprots_hgnc=exprots_hgnc[c(1:(which(rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"])[2]-1),(which(rownames(exprots_hgnc)%in%myannot$ensembl_gene_id[myannot$hgnc_id=="HGNC:30046"])[2]+1):nrow(exprots_hgnc)),]

##################CHECK BIASES########################################################
#format data for noiseq
noiseqData = readData(data = exprots_hgnc, gc = myannot[,1:2], biotype = myannot[,c(1,3)],
	factor=designExp, length=myannot[,c(1,8)])
#1)check expression bias per sample
mycountsbio = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_count_distribution_global.pdf", width = 15, height = 7)
explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:ncol(expr),200),Colv=NA)
dev.off()
#pdf("noiseqPlot_count_distribution_protein_coding.pdf", width = 15, height = 7)
#explo.plot(mycountsbio, plottype = "boxplot", samples = sample(1:889,200),toplot="protein_coding")
#dev.off()
## Biodetection plot
#mybiodetection <- dat(noiseqData, type="biodetection", factor=NULL)
#pdf("mybiodetection.pdf", width = 15, height = 7)
#explo.plot(mybiodetection)
#dev.off()
#saturation plot
mysaturation <- dat(noiseqData, k = 0, ndepth = 7, type = "saturation")
pdf("mysaturation", width = 15, height = 7)
explo.plot(mysaturation, toplot="protein_coding", 
    samples = c(1,12), yleftlim = NULL, yrightlim = NULL)
dev.off()

#2)check for low count genes
myCounts = dat(noiseqData, type = "countsbio", factor = NULL)
pdf("noiseqPlot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = sample(1:ncol(expr),200))
dev.off()
pdf("lowCountThres.pdf")
hist(rowMeans(cpm(exprots_hgnc,log=T)),ylab="genes",xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias
#each sample s is compared to the reference sample r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). Confidence intervals for the M median is
#computed by bootstrapping. If the median of M values for each comparison is not in the CI
# it means that the deviation of the sample is statistically significant. Therefore, a
# normalization procedure should be used 
mycd = dat(noiseqData, type = "cd", norm = FALSE) #slooooow
#[1] "Warning: 739 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#1081    133 
#explo.plot(mycd,samples=sample(1:ncol(expr),10))

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", factor = "definition")
pdf("GCbias.pdf")
  explo.plot(myGCcontent, samples = NULL)
 dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "definition")
pdf("lengthbias.pdf")
  explo.plot(mylenBias, samples =NULL)
dev.off()
#BUT, since the gene has the same length in all your samples, there is no need to divide by the gene length

#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = FALSE, logtransf = FALSE)
pdf("noiseqPlot_PCA_before_normalization.pdf", width = 5*2, height = 5)
par(mfrow = c(1,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "definition")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "definition")
dev.off()

#################SOLVE BIASES######################################################
#1) filter low count genes with the threshold found in BIASES(2): The CPM are used because the size 
#of the library per sample affects the power of the experiment. CPM=(counts/fragments
# sequenced)*one million. Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
countMatrixFiltered = exprots_hgnc[rowMeans(exprots_hgnc)>10,]#como Diana
#13571 features 

#2)normaliza RNA composition and then GC content coz the 1st needs raw counts
#and gives you scaling factors that can be exported to the 2nd
#factorsTMM=calcNormFactors(countMatrixFiltered)
#normTMM=cqn(countMatrixFiltered,lengths=myannot$lenght[
#  myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  sizeFactors=factorsTMM,verbose=T)
#normalisedTMMatrix=normTMM$y+normTMM$offset
#noiseqData = readData(data = normalisedTMMatrix, factors=designExp, gc = myannot[,c(1,2)])
#mycdTMM = dat(noiseqData, type = "cd", norm = T)
#"Diagnostic test: FAILED. Normalization is required to correct this bias."
#factorsUQUA=calcNormFactors(countMatrixFiltered, method="upperquartile")
#normUQUA=cqn(countMatrixFiltered,lengths=myannot$lenght[
#  myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  sizeFactors=factorsUQUA,verbose=T)
#normalisedUQUAatrix=normUQUA$y+normUQUA$offset
#noiseqData = readData(data = normalisedUQUAatrix, factors=designExp, gc = myannot[,c(1,2)])
#mycdUQUA = dat(noiseqData, type = "cd", norm = T)
#table(mycdUQUA@dat$DiagnosticTest[,  "Diagnostic Test"])
#rownames(designExp)=designExp$cases
#myDESEQ=DESeqDataSetFromMatrix(countData=countMatrixFiltered,colData=designExp,design=~definition)
#deseqFactors=estimateSizeFactors(myDESEQ)
#normDESEQ=cqn(countMatrixFiltered,lengths=myannot$lenght[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  x=myannot$percentage_gene_gc_content[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered)],
#  sizeFactors=sizeFactors(deseqFactors),verbose=T)
#myDESEQ=normDESEQ$y+normDESEQ$offset
#noiseqData = readData(data = myDESEQ, factors=designExp,gc = myannot[,c(1,2)])
#mycdDESEQ = dat(noiseqData, type = "cd", norm = T)
#table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#TMM.TP=normalisedTMMatrix[,colnames(normalisedTMMatrix)%in%designExp$barcode[designExp$shortLetterCode=="TP"]]
#TMM.NT=normalisedTMMatrix[,colnames(normalisedTMMatrix)%in%designExp$barcode[designExp$shortLetterCode=="NT"]]
#pdf("sampleBias.pdf")
#temp=sample(1:ncol(expr),10)
#explo.plot(mycd,samples=temp)
#explo.plot(mycdTMM,samples=temp)
#dev.off()


#####sigo el pipeline de Diana para normalizar porque cqn aplana demasiado la vari
#nombro las columnas de myannot y design de acuerdo a la sig función
mean10=list(M=countFilt,Annot=annot,Targets=tissue)
mydataM10EDA <- EDASeq::newSeqExpressionSet(counts=as.matrix(countFilt),featureData=annot,phenoData=tissue)
PLOTSNORMDIR ="plots"
dir.create(PLOTSNORMDIR)

getNOISeqResults <- function(step1, step2, step3, n.counts, m10.data) {
### Check the NOISEq results 
mydata <- NOISeq::readData(
  data = n.counts, 
  length = m10.data$Annot[, c("EnsemblID", "Length")], 
  biotype = m10.data$Annot[, c("EnsemblID", "Type")], 
  #     chromosome = m10.data$Annot[, c("Chr", "Start", "End")], 
  factors = m10.data$Targets[, "Group",drop=FALSE], 
  gc = m10.data$Annot[, c("EnsemblID", "GC")])
nsamples <- dim(m10.data$Targets)[1]

### Length bias 
mylengthbias <- dat(mydata, factor="Group", norm = TRUE, type="lengthbias")
l.stats.1 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[1])
l.stats.2 <- getRegressionStatistics(mylengthbias@dat$RegressionModels[2])
## GC Bias
mygcbias <- dat(mydata, factor = "Group", norm = TRUE, type ="GCbias")
gc.stats.1 <- getRegressionStatistics(mygcbias@dat$RegressionModels[1])
gc.stats.2 <- getRegressionStatistics(mygcbias@dat$RegressionModels[2])
   
#RNA Composition
myrnacomp <- dat(mydata, norm = TRUE, type="cd")
dtable <- table(myrnacomp@dat$DiagnosticTest[,  "Diagnostic Test"])
if (is.na(dtable["PASSED"])) dtable <- data.frame(PASSED = 0)

pngPlots <- c(paste(PLOTSNORMDIR, "Lenghtbias.png", sep="/"), 
              paste(PLOTSNORMDIR, "GCbias.png", sep="/"), 
              paste(PLOTSNORMDIR, "RNAComposition.png", sep="/"))

png(pngPlots[1], width=w/2, height=h/2)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()

png(pngPlots[2], width=w/2, height=h/2)
explo.plot(mygcbias, samples = NULL, toplot = "global")
dev.off()

png(pngPlots[3],width=w/2, height=h/2)
explo.plot(myrnacomp, samples = 1:12)
dev.off()

thePlots <- lapply (pngPlots, function(pngFile) {
  rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
})

plotname <- paste(step1, step2, step3, sep = "_")
png(paste(PLOTSNORMDIR, paste(plotname, "png", sep ="."), sep="/"),  height=h/2, width=w*(3/2))
par(oma = c(0, 0, 1.5, 0))
plot.new()
do.call(grid.arrange, c(thePlots,  ncol = 3, nrow=1))
mtext(plotname, outer = TRUE, cex = 1.5)
dev.off()

unlink(pngPlots)

norm.set.results <- data.frame(step1, step2, step3, 
                               l.stats.1$r2, l.stats.1$p, l.stats.2$r2, l.stats.2$p,
                               gc.stats.1$r2, gc.stats.1$p, gc.stats.2$r2, gc.stats.2$p, 
                               dtable["PASSED"], dtable["PASSED"]/nsamples)
colnames(norm.set.results) <- c("Step1", "Step2", "Step3", 
                                paste("Lenght", l.stats.1$name, "R2", sep = "."), paste("Lenght", l.stats.1$name, "p-value", sep = "."),  
                                paste("Lenght", l.stats.2$name, "R2", sep = "."), paste("Lenght", l.stats.2$name, "p-value", sep = "."), 
                                paste("GC", gc.stats.1$name, "R2", sep = "."), paste("GC", gc.stats.1$name, "p-value", sep = "."), 
                                paste("GC", gc.stats.2$name, "R2", sep = "."), paste("GC", gc.stats.2$name, "p-value", sep = "."),
                                "RNA.PassedSamples", "RNA.PassedProportion")
return(norm.set.results)
  } 
 getRegressionStatistics <- function(regressionmodel) {
    name <- names(regressionmodel)
    print(name)
    rsquared <- summary(regressionmodel[[1]])$r.squared
    print(rsquared)
    fstatistic <- summary(regressionmodel[[1]])$fstatistic
    pvalue <- signif(pf(q = fstatistic[1], df1 = fstatistic[2], df2 = fstatistic[3], lower.tail = FALSE), 2)
    return(list("name" = name, "r2" = rsquared, "p" = pvalue))
  }

  ## We try with GC normalization first
  for (gcn in gc.norm) {
    gcn.data <- withinLaneNormalization(mydataM10EDA, "GC", which = gcn)
    for (ln in lenght.norm) {
      ln.data <- withinLaneNormalization(gcn.data, "Length", which = ln)
      for (bn in between.nom) {
        if (bn == "tmm") {
          between.data <- tmm(normCounts(ln.data), long = 1000, lc = 0, k = 0)
        } else {
          between.data <- betweenLaneNormalization(ln.data, which = bn, offset = FALSE)
        }
        cat("Testing with GC normalization: ", gcn, ",  length normalization: ", ln, " and between lane normalization: ", bn, "\n")
        norm.noiseq.results <- getNOISeqResults(paste("GC", gcn, sep = "."), paste("Length", ln, sep = "."), paste("Between", bn, sep =  "."), 
                                                between.data, mean10)
        normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
      }
    }
  }
  
  ## We try with length normalization now
  for (ln in lenght.norm) {
    ln.data <- withinLaneNormalization(mydataM10EDA, "Length", which = ln)
    for (gcn in gc.norm) {
      gcn.data <- withinLaneNormalization(ln.data, "GC", which = gcn)
      for (bn in between.nom) {
        if (bn == "tmm") {
          between.data <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
        } else {
          between.data <- betweenLaneNormalization(gcn.data, which = bn, offset = FALSE)
        }
        cat("Testing with length normalization: ", ln, ", GC normalization: ", gcn, " and between lanes normalization: ", bn, "\n")
        norm.noiseq.results <- getNOISeqResults(paste("Length", ln, sep = "."), paste("GC", gcn, sep = "."), paste("Between", bn, sep =  "."),
                                                counts(between.data), mean10)
        normalization.results <- rbind(normalization.results, norm.noiseq.results)                      
      }
    }
  }

 
#######################starting classification########################################################
subannot=getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene"),values=rownames(TMM.TP), mart=mart)
colnames(subannot)=c("probe","NCBI.gene.symbol","EntrezGene.ID")#this colnames are needed
subannot=subannot[subannot$probe%in%rownames(TMM.TP),]
subannot=subannot[!duplicated(subannot$probe),]
subannot=subannot[order(match(subannot$probe,rownames(TMM.TP))),]

#49/50 probes are used for clustering
subtypes=molecular.subtyping(sbt.model="pam50",data=t(log2(TMM.TP)),annot=subannot,do.mapping=T)
subtipos=TCGA_MolecularSubtype(colnames(TMM.TP))
table(subtipos$gender)
#female   male 
# 1201     13 
subtipos=subtipos$subtypes[,c(4,2)]
temp=cbind(colnames(TMM.TP)[!colnames(TMM.TP)%in%subtipos$barcodes],NA)
colnames(temp)=colnames(subtipos)
subtipos=rbind(subtipos,temp)
subtipos=cbind(subtipos,subtypes$subtype[order(match(names(subtypes$subtype),subtipos$barcode))])
table(subtipos[,2:3],useNA="ifany") 
#                  X
#subtype_PAM50.mRNA Basal Her2 LumB LumA Normal
#     Basal-like       97    0    2    2      0    0
#     HER2-enriched     1   52    5    0      0    0
#     Luminal A         2    6   42  180      2    0
#     Luminal B         1    8  112    4      0    0
#     Normal-like       4    0    2    0      2    0
#     <NA>            100   55  152  251     20    0

#############################solve batch effect#######################################################
temp=cbind(designExp$barcode[!designExp$barcode%in%subtipos[,1]],"normal","normal","normal","normal")
subtipos=rbind(subtipos,temp)

noiseqData = readData(data = normalisedTMMatrix, factors=subtipos, gc = myannot[,c(1,2)])
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = FALSE)
TMMARSyn=ARSyNseq(noiseqData, factor = "genefu", batch = F, norm = "n",  logtransf = FALSE)
myPCA1 = dat(TMM.ARSyn, type = "PCA", norm = T)
pdf("pre-postArsynPCA.pdf")
par(mfrow=c(2,2))
explo.plot(myPCA, samples = c(1,2), plottype = "scores", factor = "genefu")
explo.plot(myPCA, samples = c(1,3), plottype = "scores", factor = "genefu")
explo.plot(myPCA1, samples = c(1,2), plottype = "scores", factor = "genefu")
explo.plot(myPCA1, samples = c(1,3), plottype = "scores", factor = "genefu")
dev.off()

TMMArsyn.NT=exprs(TMMARSyn)[,colnames(exprs(TMMARSyn))%in%subtipos[subtipos[,3]=="normal",1]]
TMMArsyn.TP=exprs(TMMARSyn)[,colnames(exprs(TMMARSyn))%in%subtipos[subtipos[,3]!="normal",1]]

##############final classification based on genefu's alredy classified samples##########################
rownames(subannot)=subannot$probe
obje<-PAM50(exprs=log2(TMMArsyn.TP),annotation=subannot)
obje=filtrate(obje,verbose=T)
obje=classify(obje,verbose=T,std="scale")
obje=permutate(obje,verbose=T,keep=T,nPerm=10000, pCutoff=0.01, where="fdr",corCutoff=0.1,seed=1234567890,
  BPPARAM=MulticoreParam(progressbar=TRUE))
#49/50 probes are used for clustering
temp=permutation(obje)$subtype
subtipos[subtipos[,3]!="normal",5]=temp$Subtype[order(match(rownames(temp),subtipos[,1]))]
#table(subtipos[,c(2,5)],useNA="always")
#               pbcmc2
#tcga            Ambiguous Basal Her2 LumA LumB Normal <NA>
#  Basal-like            0    97    0    2    1      0    1
#  HER2-enriched         0     1   52    0    5      0    0
#  Luminal A             1     2    6  180   37      2    4
#  Luminal B             2     1    8    4   93      0   17
#  Normal-like           0     4    0    0    2      2    0
#  <NA>                  6   100   54  250  131     19   18


#final transcript composition bias
noiseqData = readData(data = exprs(TMM.ARSyn.TP), factors=subtipos, gc = myannot[,c(1,2)])
mycdnoBatch=dat(noiseqData,type="cd",norm=T)
table(mycdnoBatch@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   496    605 

#final expression dataset
subtipos=subtipos[order(match(subtipos[,1],designExp$barcode)),]
subtipos=cbind(designExp,subtipos[,2:5])
#sólo me quedo con una columna por paciente, como están normalizados, según yo no me importa cual agarre
TMMArsyn=cbind(TMMArsyn.TP[,!duplicated(substr(colnames(TMMArsyn.TP),1,12))],TMMArsyn.NT)
TMMArsyn=lapply(c("Basal","Her2","LumA","LumB","normal"),function(x) TMMArsyn[,colnames(TMMArsyn)%in%subtipos$barcode[subtipos$pbcmc2==x]])
names(TMMArsyn)=c("Basal","Her2","LumA","LumB","normal")
save(TMMArsyn,subtipos,file="subtiTMMArsyn.RData")

############################--------------------FIN--------------##########################

#check biotype bias per tissue
#mybiodetection <- dat(noiseqData, k = 0, type = "biodetection", factor = "tissue.definition")
#pdf("ProteinCodingComparison.pdf",height=7,width=7)
#par(oma=c(5,1,1,1))
#explo.plot(mybiodetection,samples=c(1,2),toplot="protein_coding",plottype="comparison")
#[1] "Percentage of protein_coding biotype in each sample:"
#Primary solid Tumor Solid Tissue Normal 
#            35.4141             35.5303 
#[1] "Confidence interval at 95% for the difference of percentages: Primary solid Tumor - Solid Tissue Normal"
#[1] -0.6746  0.4421
#[1] "The percentage of this biotype is NOT significantly different for these two samples (p-value = 0.6868 )."


