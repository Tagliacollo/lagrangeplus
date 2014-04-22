#!/bin/sh -x
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -N lagrangeML.sh
#PBS -o stdout.txt
#PBS -e stderr.txt
#PBS -A loni_sylv2012
#PBS -S /bin/sh

                 MODEL=$PBS_O_WORKDIR"/lagrangeML.py"
            configFile='__FAMILY__.config'
outputFilenameTemplate='__FAMILY_____INDEX_____REP___results.csv'

count=0
for node in `cat $PBS_NODEFILE`; do
    outputFilename=`echo $outputFilenameTemplate | sed s/__REP__/$count/g`
    (ssh $node $MODEL -w $PBS_O_WORKDIR -m $configFile -o $outputFilename )&
    count=$(( $count + 1 ))
done

echo "Waiting for child processes to finish"
wait
echo "Done waiting"

