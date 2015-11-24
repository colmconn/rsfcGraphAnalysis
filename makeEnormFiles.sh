#!/bin/bash

## set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts


# taskFile=taskFile-enorm.txt
# cat /dev/null > $taskFile

# subjects=$( cd $DATA; ls -1d *_A )
# for ss in $subjects ; do
#     ddir=/data/sanDiego/rsfcGraphAnalysis/data/$ss/rsfcPreprocessed/tmp
#     if [[ -f $ddir/${ss}funcon.zp.float.despike_tsh_vr_motion.1D ]] ; then
# 	echo "cd $ddir; 1d_tool.py -infile ${ss}funcon.zp.float.despike_tsh_vr_motion.1D -set_nruns 1 -collapse_cols euclidean_norm -write ${ss}.e.norm.1D" >> $taskFile
#     fi
# done
		 
# parallel  --wd $( pwd ) -j50 < $taskFile

subjects=$( cat $GROUP_DATA/subjectOrder.mddAndCtrl.L_BLA.weight.3mm.csv | sed 1d )

averageEnormFile=$GROUP_DATA/restingstate.mddAndCtrl.motion.enorm.csv
echo "subject, enorm" > $averageEnormFile
for ss in $subjects ; do
    ddir=/data/sanDiego/rsfcGraphAnalysis/data/$ss/rsfcPreprocessed/tmp
    enormFile=$ddir/${ss}.e.norm.1D
    if [[ -f $enormFile ]] ; then
	enorm=$( cat $enormFile | awk '{s+=$1}END{print s/NR}' )
    else
	enorm="NA"
    fi
    echo "$ss, $enorm" >> $averageEnormFile
done
