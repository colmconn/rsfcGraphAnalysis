#!/bin/bash

#set -x

## subjects=$( cat ../data/config/subject_list_from_matthew.txt )

SUBJECTS_DIR=/data/sanDiego/freesurferAnalysis/data

subjects=$( cd ../data ; ls -d *_A* )

#subjects="105_A"
echo $subjects

problemBrains=""
if [[ ! -d ../data/vbmFromFreesurfer.sobp2017/struc ]] ;then
    mkdir -p ../data/vbmFromFreesurfer.sobp2017/struc

fi

## the following is the list of medicated subjects that Laura put together for the T1 + telomere paper.
# it is taken from /data/sanDiego/structuralTelomereAnalysis/scripts/analyseHippocampalVolumes_Laura.r
## Now remove subjects that are NOT medication-naive (Prozac, Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbatta)
#mgd <- mgd[! mgd$subject %in% c(130, 132, 318, 319, 320, 322, 323, 324, 325, 329, 345, 376, 384), ]

## Now remove subjects that were on ADHD medication (Focolin, Stratera)
#mgd <- mgd[! mgd$subject %in% c(333, 349), ]


for subject in ${subjects} ; do

    if [[ $subject != 130_A ]] && \
	   [[ $subject != 132_A ]] && \
	   [[ $subject != 318_A ]] && \
	   [[ $subject != 319_A ]] && \
	   [[ $subject != 320_A ]] && \
	   [[ $subject != 322_A ]] && \
	   [[ $subject != 323_A ]] && \
	   [[ $subject != 324_A ]] && \
	   [[ $subject != 325_A ]] && \
	   [[ $subject != 329_A ]] && \
	   [[ $subject != 333_A ]] && \
	   [[ $subject != 345_A ]] && \
	   [[ $subject != 349_A ]] && \
	   [[ $subject != 376_A ]] && \
	   [[ $subject != 384_A ]] ; then

	## original list of exclusions that includes medicated and other subjects. The reason for the exclusionof the other subjects is unclear 
	# if [[ $subject != 119_A ]] && \
	    # 	   [[ $subject != 346_A ]] && \
	    # 	   [[ $subject != 370_A ]] && \
	    # 	   [[ $subject != 130_A ]] && \
	    # 	   [[ $subject != 132_A ]] && \
	    # 	   [[ $subject != 319_A ]] && \
	    # 	   [[ $subject != 320_A ]] && \
	    # 	   [[ $subject != 322_A ]] && \
	    # 	   [[ $subject != 323_A ]] && \
	    # 	   [[ $subject != 325_A ]] && \
	    # 	   [[ $subject != 329_A ]] && \
	    # 	   [[ $subject != 333_A ]] && \
	    # 	   [[ $subject != 376_A ]] && \
	    # 	   [[ $subject != 311_A ]] && \
	    # 	   [[ $subject != 378_A ]] ; then 
	
	if [[ ! -f $SUBJECTS_DIR/${subject}/mri/brainmask.mgz ]] ; then
	    problemBrains="$subject $problemBrains"
	else
	    # if [[ ! -f ../${subject}/anat/${subject}.anat_struc.std.nii.gz ]] ;then
	    # 	echo "No such file ../${subject}/anat/${subject}.anat_struc.std.nii.gz"
	    # 	exit
	    # fi
	    ( cd ../data/vbmFromFreesurfer.sobp2017      ;  ln -sf ../${subject}/anat/${subject}.anat_struc.std.nii.gz   ${subject}.anat.nii.gz )
	    ( cd ../data/vbmFromFreesurfer.sobp2017/struc; \
	      mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/T1.mgz ${subject}.anat_struc.RSP.nii; \
	      3dresample -orient RPI -inset ${subject}.anat_struc.RSP.nii -prefix ${subject}.anat_struc.nii ; \
	      rm -f ${subject}.anat_struc.RSP.nii )
	    ( cd ../data/vbmFromFreesurfer.sobp2017/struc; \
	      mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/brainmask.mgz  ${subject}.anat_struc_brain.RSP.nii ; \
	      3dresample -orient RPI -inset ${subject}.anat_struc_brain.RSP.nii -prefix ${subject}.anat_struc_brain.nii ; \
	      rm -f ${subject}.anat_struc_brain.RSP.nii )
	fi
    else
	echo "Skipping $subject"
    fi
done
		    
echo "Print problem brains: $problemBrains"
