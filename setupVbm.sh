#!/bin/bash

# set -x

## subjects=$( cat ../data/config/subject_list_from_matthew.txt )

ctrlSubjects=$( cat ../data/vbm.subject.list.from.matthew.n111/rs114a.ncl.txt )
mddSubjects=$( cat ../data/vbm.subject.list.from.matthew.n111/rs114a.mdd.txt )

## subjects="${ctrlSubjects} ${mddSubjects}"


## subjects=$( cd ../data ; ls -d *_A* )
echo $subjects

problemBrains=""
if [[ ! -d ../data/vbm.subject.list.from.matthew.n111/struc ]] ;then
    mkdir -p ../data/vbm.subject.list.from.matthew.n111/struc

fi

for subject in ${ctrlSubjects} ; do

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
	
	if [[ ! -f ../data/${subject}/anat/${subject}.anat_struc.std.nii.gz ]] ; then
	    problemBrains="$subject $problemBrains"
	else
	    ( cd ../data/vbm.subject.list.from.matthew.n111      ;  ln -sf ../${subject}/anat/${subject}.anat_struc.std.nii.gz   NCL.${subject}.anat.nii.gz )
	    ( cd ../data/vbm.subject.list.from.matthew.n111/struc;  ln -sf ../../${subject}/anat/${subject}.anat_struc.std.nii.gz   NCL.${subject}.anat_struc.nii.gz )
	    ( cd ../data/vbm.subject.list.from.matthew.n111/struc;  ln -sf ../../${subject}/anat/${subject}.anat_struc_brain.std.nii.gz  NCL.${subject}.anat_struc_brain.nii.gz )
	fi
    # else
    # 	echo "Skipping $subject"
    # fi
done

      
echo "Print problem brains: $problemBrains"

for subject in ${mddSubjects} ; do

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
	
	if [[ ! -f ../data/${subject}/anat/${subject}.anat_struc.std.nii.gz ]] ; then
	    problemBrains="$subject $problemBrains"
	else
	    ( cd ../data/vbm.subject.list.from.matthew.n111      ;  ln -sf ../${subject}/anat/${subject}.anat_struc.std.nii.gz   MDD.${subject}.anat.nii.gz )
	    ( cd ../data/vbm.subject.list.from.matthew.n111/struc;  ln -sf ../../${subject}/anat/${subject}.anat_struc.std.nii.gz   MDD.${subject}.anat_struc.nii.gz )
	    ( cd ../data/vbm.subject.list.from.matthew.n111/struc;  ln -sf ../../${subject}/anat/${subject}.anat_struc_brain.std.nii.gz  MDD.${subject}.anat_struc_brain.nii.gz )
	fi
    # else
    # 	echo "Skipping $subject"
    # fi
done

      
echo "Print problem brains: $problemBrains"
