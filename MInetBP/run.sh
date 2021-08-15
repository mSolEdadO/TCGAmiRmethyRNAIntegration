#!/usr/bin/env bash

BP=$1
echo "Building data for $BP"
cd data/
Rscript ../BP.R $BP
if [[ ! -f "$BP" ]];then exit 1; fi
noncpgs=$(echo "$BP.txt")
cat $BP>$noncpgs
grep -v 'probes' Her2.TFs >> $noncpgs
grep -v 'probes' Her2.miRNAs >> $noncpgs
cpgs=$(echo "$BP-cpgs.txt")
cat $BP>$cpgs
grep -v 'probes' Her2.CpGs>>$cpgs

echo "Building submit file"
cd ../
Rscript submit_file.R $BP
sub=$(echo "$BP.sub")
cat $sub|perl -pe 's/\t/\n/g'>temp
cat condor.header temp>$sub

condor_submit $sub
echo "Waiting for aracne to end" 
aim=$(wc -l bp.txt|cut -d' ' -f1)
echo "$aim" #107
aim=$(expr $aim + $aim )
echo "$aim" #214
done=$(ls data/|grep -c 'adj')
echo "$done" #0
while [ $done -le $aim ]; do 
	done=$(ls data/|grep -c 'adj');
	echo "$done";
	sleep 5;done #se queda pasmado, se produce adj aunque no halla interacciones? quiza estas esperando algo que no va a pasar
echo "Building sif file"	
line=$(echo "bin/adj2sif data/$BP_*.adj>$BP.sif")
eval $line

echo "Removing intemediate files"
#rm *output *log *error $sub
#rm data/*adj data/*txt data/GO\:*
