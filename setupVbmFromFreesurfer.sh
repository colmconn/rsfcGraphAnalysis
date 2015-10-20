#!/bin/bash

#set -x

## subjects=$( cat ../data/config/subject_list_from_matthew.txt )

SUBJECTS_DIR=/data/sanDiego/freesurferAnalysis/data

subjects=$( cd ../data ; ls -d *_A* )

#subjects="105_A"
echo $subjects

problemBrains=""
if [[ ! -d ../data/vbmFromFreesurfer/struc ]] ;then
    mkdir -p ../data/vbmFromFreesurfer/struc

fi

for subject in ${subjects} ; do

    if [[ $subject != 119_A ]] && \
	   [[ $subject != 346_A ]] && \
	   [[ $subject != 370_A ]] && \
	   [[ $subject != 130_A ]] && \
	   [[ $subject != 132_A ]] && \
	   [[ $subject != 319_A ]] && \
	   [[ $subject != 320_A ]] && \
	   [[ $subject != 322_A ]] && \
	   [[ $subject != 323_A ]] && \
	   [[ $subject != 325_A ]] && \
	   [[ $subject != 329_A ]] && \
	   [[ $subject != 333_A ]] && \
	   [[ $subject != 376_A ]] && \
	   [[ $subject != 311_A ]] && \
	   [[ $subject != 378_A ]] ; then 
	
	if [[ ! -f $SUBJECTS_DIR/${subject}/mri/brainmask.mgz ]] ; then
	    problemBrains="$subject $problemBrains"
	else
	    ( cd ../data/vbmFromFreesurfer      ;  ln -sf ../${subject}/anat/${subject}.anat_struc.std.nii.gz   ${subject}.anat.nii.gz )
	    ( cd ../data/vbmFromFreesurfer/struc; \
	      mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/T1.mgz ${subject}.anat_struc.RSP.nii; \
	      3dresample -orient RPI -inset ${subject}.anat_struc.RSP.nii -prefix ${subject}.anat_struc.nii ; \
	      rm -f ${subject}.anat_struc.RSP.nii )
	    ( cd ../data/vbmFromFreesurfer/struc; \
	      mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/brainmask.mgz  ${subject}.anat_struc_brain.RSP.nii ; \
	      3dresample -orient RPI -inset ${subject}.anat_struc_brain.RSP.nii -prefix ${subject}.anat_struc_brain.nii ; \
	      rm -f ${subject}.anat_struc_brain.RSP.nii )
	fi
    else
	echo "Skipping $subject"
    fi
done
		    
echo "Print problem brains: $problemBrains"
