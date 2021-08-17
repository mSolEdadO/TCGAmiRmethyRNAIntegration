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
if [[ ! -f "$sub" ]];then exit 1; fi
cat $sub|perl -pe 's/\t/\n/g'>temp
cat condor.header temp>$sub

condor_submit $sub
echo "Waiting for aracne to end" 
aim=$(condor_q|grep 'jobs'|cut -d' ' -f1)
while [ $aim -gt 0 ]; do 
	echo "$aim";
	aim=$(condor_q|grep 'jobs'|cut -d' ' -f1);
	sleep 5;
done
echo "Building sif file"	
line=$(echo "bin/adj2sif data/$BP*.adj>$BP.sif")
eval $line
echo "Removing intemediate files"
rm *output *log *error $sub bp.txt temp
rm data/GO\:*
####does not come back after ending
