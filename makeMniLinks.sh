#!/bin/bash

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data

function makeMniLink {
    local directory=$1

    echo "*** Making MNI links in $directory"
    
    owd=$( pwd ) 
    cd $directory

    if [[ ! -f MNI152_T1_1mm.nii.gz  ]] ; then
	ln -sf $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz
    fi

    if [[ ! -f MNI152_T1_1mm_brain.nii.gz ]] ; then
	ln -sf $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz
    fi

    symlinks -c ./

    cd $owd
}

groupResultsDirs="$( ls -d $DATA/Group.results* )"

for dd in ${groupResultsDirs} ; do
    makeMniLink $dd
done

# makeMniLink $DATA/Group.results
# makeMniLink $DATA/Group.results.AAL

# regressionVariables="CDRSR.diff MASC.tscore.diff CGAS.diff RADS.Total.Tscore.diff"

# GROUP_RESULTS=$DATA/Group.results.$regressionVariable.withAandC

# for dd in $regressionVariables ; do
#     makeMniLink $DATA/Group.results.$dd.withAandC
# done
    
    
