#!/bin/bash

# set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data

subjects="$( cat $GROUP_DATA/subjectOrder.mddAndCtrl.HO.L.amygdala.3mm.csv | sed 1d )"

[[ ! -d $DATA/ica ]] && mkdir $DATA/ica
cd $DATA/ica

cat /dev/null > anatFilesList.txt
cat /dev/null > restingFilesList.txt
cat /dev/null > parallel-Taskfile.txt

echo "Subject,headFile,anatFile,restingFile"
for subject in $subjects ; do

    headFile=$DATA/$subject/anat/$subject.anat_struc.std.nii.gz
    anatFile=$DATA/$subject/anat/$subject.anat_struc_brain.std.nii.gz    
    restingFile=$DATA/$subject/${subject}funcon+orig.HEAD

    if [[ -f $headFile ]] && [[ -f $anatFile ]] && [[ $restingFile ]] ; then
     	 #( cd $DATA/ica; ln -f $headFile $( echo ${headFile##*/} | sed "s/.std//" | sed "s/_A./_A_/" ) )
	 #( cd $DATA/ica; ln -f $anatFile $( echo ${anatFile##*/} | sed "s/.std//" | sed "s/_A./_A_/" ) )
	# 	echo "( cd $DATA/ica; 3dresample -orient RPI -inset $restingFile -prefix ${subject}_restingstate.nii )" >> parallel-Taskfile.txt
	
	# 	#( cd $DATA/ica; 3dcopy $restingFile $subject.restingstate.nii )

     	echo "$DATA/ica/$( echo ${anatFile##*/} | sed "s/.std//" | sed "s/_A./_A_/" )" >> anatFilesList.txt
	echo "$DATA/ica/${subject}_restingstate.nii.gz" >> restingFilesList.txt
	
	# 	echo "*** Including $subject"
	# else 
	# 	echo "*** Skipping $subject"
    fi
    
    csvLine="$subject"
    if [[ -f "$( echo ${headFile##*/} | sed "s/.std//" | sed "s/_A./_A_/" )" ]] ; then
	#echo "*** $subject head file good"
	csvLine="$csvLine,good"
    else
	csvLine="$csvLine,bad"	
    fi
    
    if [[ -f "$( echo ${anatFile##*/} | sed "s/.std//" | sed "s/_A./_A_/" )"  ]] ; then
	#echo "*** $subject brain anat file good"
	csvLine="$csvLine,good"
    else
	csvLine="$csvLine,bad"		
    fi
    
    if [[ -f ${subject}_restingstate.nii.gz ]] ; then
	#echo "*** $subject resting state file good"
	csvLine="$csvLine,good"
    else
	csvLine="$csvLine,bad"		
    fi

    echo $csvLine

done

symlinks -c ./
# parallel --wd $( pwd ) -j16 < parallel-Taskfile.txt
