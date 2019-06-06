#load files
files=list.files()
#files=files[grep("condor",files,invert=T)]
files=files[5:length(files)]
subtis=unique(sapply(strsplit(files,'.',fixed=T),function(x) x[2]))
subtis=lapply(subtis,function(x) files[grep(x,files)])
names(subtis)=c("Basal","Her2","LumB","normal","LumA")
modelos=lapply(subtis,function(x) lapply(x,function(y) read.table(y,header=T)))

#check fitted parameters
params=lapply(modelos,function(x) sapply(x,colnames))
params=lapply(params,function(x) gsub('X','',x))
params=lapply(params,function(x) t(do.call(cbind,strsplit(x,'_'))))
png("alpha.png")
 par(mar=c(2,10,2,2))
 barplot(t(table(temp[,c(1,3)])),beside=T,las=2,horiz=T,xlim=c(0,30))
dev.off()
png("lambda.png")
 par(mar=c(2,10,2,2))
 barplot(t(table(temp[,2:3])),beside=T,las=2,horiz=T,xlim=c(0,25))
dev.off()

#sif per subtype
subtis=lapply(subtis,function(x) sapply(strsplit(x,'.',fixed=T),function(y) y[1]))
modelos=lapply(modelos,function(y) lapply(y,function(x) as.matrix(x)[which(x>0),]))
modelos=lapply(1:5,function(x) sapply(1:length(modelos[[x]]),function(y) 
    cbind(subtis[[x]][y],names(modelos[[x]][[y]]),modelos[[x]][[y]])))
modelos=lapply(modelos,function(x) do.call(rbind,x))
sapply(modelos,function(x) summary(as.numeric(x[grep("(Intercept)",x[,2]),3])))
            Basal      Her2      LumB    normal      LumA
Min.    8.615e-18 1.052e-16 4.303e-17 6.979e-17 1.052e-16
1st Qu. 6.219e-16 1.109e-15 3.279e-16 5.092e-16 1.109e-15
Median  1.125e-15 1.527e-15 6.452e-16 1.145e-15 1.527e-15
Mean    1.597e-15 1.474e-15 1.133e-15 1.407e-15 1.474e-15
3rd Qu. 2.473e-15 2.144e-15 1.934e-15 1.947e-15 2.144e-15
Max.    4.190e-15 2.860e-15 3.105e-15 3.668e-15 2.860e-15
sapply(modelos,function(x) summary(as.numeric(x[grep("(Intercept)",x[,2],invert=T),3])))
            Basal      Her2      LumB    normal      LumA
Min.    6.881e-10 2.636e-08 3.198e-07 4.061e-11 2.636e-08
1st Qu. 5.538e-05 3.677e-04 7.383e-04 6.082e-05 3.677e-04
Median  1.373e-04 1.009e-03 2.297e-03 1.509e-04 1.009e-03
Mean    9.030e-04 4.510e-03 8.468e-03 9.555e-04 4.510e-03
3rd Qu. 2.901e-04 2.752e-03 6.914e-03 3.191e-04 2.752e-03
Max.    8.175e-01 7.335e-01 7.671e-01 3.909e-01 7.335e-01
#modelos=sapply(modelos,function(x) x[grep("(Intercept)",x[,2],invert=T),])

#get node chr, hgnc_id...