#!/bin/bash

## set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$RESULTS/Group.results
MDD_STANDARD=$ROOT/standard
scriptDir=$ROOT/scripts


subjects="$( sed 1d < $GROUP_DATA/subjectOrder.mddAndCtrl.L_whole_amygdala.3mm.csv )"
seedList=$DATA/config/juelich_whole_amygdala_seeds.txt
seeds=$( eval echo $( cat $seedList ) )
seeds="$MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz $seeds"

## subjects=105_A

for subject in ${subjects}; do

    echo "*** Warping TSNR to MNI152 space"

    cd $DATA/$subject/rsfcPreprocessed
    if [[ ! -f ${subject}.tsnr.MNI.nii.gz ]] ; then 
	3dcopy ${subject}.pm.tsnr+orig.HEAD ${subject}.pm.tsnr.nii
	
	3dresample -master ../anat/${subject}.anat_struc.std.nii.gz -prefix ./$subject.tsnr.resampled.nii \
		   -inset ${subject}.func.native.std.nii.gz
	
	flirt -in ${subject}.pm.tsnr.nii -ref ${subject}.tsnr.resampled \
	      -out ${subject}.tsnr2anat.flirt -omat ${subject}.tsnr2anat.flirt.mat
	
	applywarp						\
	    --ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz		\
	    --in=${subject}.pm.tsnr.nii.gz			\
	    --warp=../anat/${subject}.std.2.MNI.warpcoef	\
	    --premat=${subject}.tsnr2anat.flirt.mat		\
	    --out=${subject}.tsnr.MNI
    fi
    
    for seed in $seeds ; do
	
    	seedName=${seed##*/}
    	if echo $seedName | grep -q "nii" ; then 
    	    seedName=${seedName%%.nii*}
    	else 
    	    seedName=${seedName%%+*}
    	fi
    	echo "*** Extracting TSNR for subject $subject from $seed"

	3dmaskave -quiet -mask $seed ${subject}.tsnr.MNI.nii.gz > ${subject}.${seedName}.tsnr.1D

    done						
done
