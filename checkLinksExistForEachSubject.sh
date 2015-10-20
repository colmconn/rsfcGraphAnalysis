#!/bin/bash

##set -x 

cd /data/sanDiego/rsfcGraphAnalysis/data

subjects=$( ls -d [0-9][0-9][0-9]_[AC]* ) 
##subjects="311_A"

echo "subject,anat HEAD,anat BRIK,funcon HEAD,funcon BRIK,all Links OK?, cleanEPI.MNI exists, #ventricle voxels in mask, #wm voxels in mask"
for ss in $subjects ; do
    echo -n "$ss,"
    nnFiles=0
    if [[ -L $ss/${ss}+orig.HEAD ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ -L $ss/${ss}+orig.BRIK.gz ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi


    if [[ -L $ss/${ss}funcon+orig.HEAD ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ -L $ss/${ss}funcon+orig.BRIK.gz ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ $nnFiles -eq 4 ]] ; then  
	echo -n "TRUE,"
    else 
	echo -n "FALSE,"
    fi

    if [[ -f $ss/rsfcPreprocessed/${ss}.pm.cleanEPI.MNI.nii.gz ]] ; then 
	echo -n "TRUE,"
    else
	echo -n "FALSE,"
    fi

    if [[ -f $ss/rsfcPreprocessed/tmp/mask.vent.resample.erode+orig.HEAD ]] ; then
	nvoxels="$( cd $ss/rsfcPreprocessed/tmp/ ;  3dclust  -isomerge  -quiet 0 0 mask.vent.resample.erode+orig. 2>/dev/null| awk '{print $1}' )"
	echo "$nvoxels" | grep -q "#" && echo -n "No voxels," || echo -n "$nvoxels, "
    else
	echo -n "No File,"
    fi

    if [[ -f $ss/rsfcPreprocessed/tmp/mask.WM.resample.erode+orig.HEAD ]] ; then
	nvoxels="$( cd $ss/rsfcPreprocessed/tmp/ ;  3dclust  -isomerge  -quiet 0 0 mask.WM.resample.erode+orig. 2>/dev/null| awk '{print $1}' )"
	echo "$nvoxels" | grep -q "#" && echo "No voxels" || echo "$nvoxels"
    else
	echo "No File"
    fi


done


#{118_C,120_C,149_C,150_C,301_C,311_A,315_C,364_A,364_C}

#118_C 120_C 149_C 150_C 301_C 311_A 315_C 364_A 364_C
