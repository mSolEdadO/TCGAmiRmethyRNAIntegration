library(data.table)

#filter out links with BPs
i=sapply(regus,function(x) names(x)[1])
i=sapply(1:5,function(x) which(BPenriched[[x]]$V1==i[x])[1]-1)
temp=lapply(1:5,function(x) BPenriched[[x]][1:i[x],])
#get the index of each kind of interaction
i=lapply(temp,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
j=unique(i[[1]])
#min MI on BPenriched nets
minMI=sapply(1:5,function(x) sapply(j,function(y) 
	round(min(as.numeric(temp[[x]]$V2)[i[[x]]==y]),4)))
colnames(minMI)=names(top)

##read sifs at different p.value thresholds
files=list.files()
files=files[grep("sif",files)]
sif=lapply(files,fread)
names(sif)=gsub(".sif","",files)
#index per interaction type
i=lapply(sif,function(x) paste(substr(x$V1,1,1),substr(x$V3,1,1)))
temp=sapply(1:length(i),function(x) summary(sif[[x]]$V2[i[[x]]=="h E"]))
colnames(temp)=names(sif)
#which contains the MI values on nets
