library(data.table)
her2=fread("data/her2.txt")

expre=as.matrix(her2[c(1:447,385021:398925),2:76])
rownames(expre)=her2$probes[c(1:447,385021:398925)]
methy=as.matrix(her2[,2:76])
rownames(methy)=her2$probes

i=seq(1,nrow(methy),12000)
i=i[1:(length(i)-1)]
chunks=lapply(i,function(x) methy[x:(x++11999),])
chunks[[33]]=methy[(i[length(i)]+11999):nrow(methy),]
lapply(1:33,function(x) 
	write.table(chunks[[x]],paste(paste("data/sub",x,sep=""),"txt",sep='.'),sep='\t',quote=F))

#$ sed -n '12,6007334 p' condor.sub|perl -pe 's/.*?(ENSG|hsa).*?\n//g;unless(/queue/){s/\n/ñ/;s/timeñqueue//g}'>condorsub.edited
#$ for x in {1..33};do y=$(echo 'sub'$x);sed -n '1,11 p' condor.sub>chunks/$y;grep -w $y condorsub.edited|perl -pe 's/ñ/\n/g'>>chunks/$y;done
#$ for x in $(ls );do condor_submit $x;done