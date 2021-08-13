#!/usr/bin/env bash

BP=$1
echo "Building data for $BP"
Rscript BP.R $BP
noncpgs=$(echo "data/$BP.txt")
cat ${BP}>${noncpgs}
grep -v 'probes' Her2.TFs>>${noncpgs}
grep -v 'probes' Her2.miRNAs>>${noncpgs}
cpgs=$(echo "data/$BP-cpgs.txt")
cat ${BP}>${cpgs}
grep -v 'probes' Her2.CpGs>>${cpgs}

echo "Building submit file"
Rscript submit_file.R $BP
sub=$(echo "$BP.sub")
cat sub|perl -pe 's/\t/\n/g'>temp
cat condor.header temp>${sub}

echo "Submiting jobs"
condor_submit $sub

#when all jobs are done, get the sif with comments of kernel size...
#erase intemediate files
a=$(condor_q -run)
echo $a|wc
while $a>1; wait
line=$(echo "bin/adj2sif data/$BP_*.adj>$BP.sif")
eval $line
rm *output *log *error
rm data/*adj