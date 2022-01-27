methy=fread("../../../Downloads/HumanMethylation450_15017482_v1-2.csv",skip=7,sep=',')
methy=methy[,c(1,22,24,26,29,30)]
desglose=apply(methy,1,function(x) cbind(unlist(strsplit(x[2],";")),unlist(strsplit(x[3],";"))))
desglose=lapply(desglose,unique)
desglose=do.call(rbind,lapply(1:length(desglose),function(x) cbind(methy[x,1],desglose[[x]])))
colnames(desglose)[2:3]=c("gene","position")
write.table(desglose,"../ini/MapMethyProbeId.tsv",quote=F,row.names=F,sep='\t')
