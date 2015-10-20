#!/bin/bash

#set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results


function linkFiles {
    local regressionVariable=$1

    suffix=Group.results.$regressionVariable.withAandC
    
    GROUP_RESULTS=$DATA/$suffix
    if [[ ! -d $GROUP_RESULTS ]] ; then
	echo "No such directory: $GROUP_RESULTS"
	break
    fi

    cd $GROUP_RESULTS

    for ff in ../20percentMotionThresholdWholeBrainClusterCorrection/$suffix/restingstate.mddOnly.*.3mm.$regressionVariable.rlm.bucket.* ; do
	ln -sf $ff ./
    done
}


regressionVariables="CDRSR.diff MASC.tscore.diff CGAS.diff RADS.Total.Tscore.diff"

for regressionVariable in $regressionVariables ; do

    linkFiles $regressionVariable
done
    
