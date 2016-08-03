#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results.CDRS.t.score.scaled.diff
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

nctrls=$( cat $GROUP_DATA/subjectOrder.ctrlOnly.R_whole_amygdala.3mm.csv | sed 1d | wc -l )
ctrls=$( cat $GROUP_DATA/subjectOrder.ctrlOnly.R_whole_amygdala.3mm.csv | sed 1d )
seedName="R_whole_amygdala.3mm"

zscoreCount=0
zscoreList=""
for subject in $ctrls ; do

    if [[ -f $DATA/$subject/rsfc/${seedName}/${seedName}.z-score+tlrc.HEAD ]] ; then 
	zscoreList="$zscoreList $DATA/$subject/rsfc/${seedName}/${seedName}.z-score+tlrc.HEAD"
	(( zscoreCount=zscoreCount+1 ))
    fi
done

if [[ $nctrls != $zscoreCount ]] ;then
    echo "*** The number of controls $nctrls is not equal to the number of zscore files ($zscoreCount)"
    exit
fi

cd $GROUP_RESULTS
3dROIstats -mask_f2short \
	   -mask R_whole_amygdala.3mm.and.CDRS.t.score.scaled.diff.R_insula.only+tlrc.HEAD \
	   ${zscoreList} > r_amy_r_insula_ctrlsOnly.1D
