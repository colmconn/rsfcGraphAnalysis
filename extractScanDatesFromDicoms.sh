#!/bin/bash

## set -x

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

dicomRoot=../../

subjectOrderFiles="new.mdd.subjectList.with.CDRS.t.score.AandC.txt new.ncl.subjectList.with.CDRS.t.score.AandC.txt"

echo "subject,group,timepoint,date"
for subjectOrder in $subjectOrderFiles ; do
    group="$( echo $subjectOrder | cut -d. -f2 )"
    subjects=$( cat ../data/config/$subjectOrder )

    for subject in $subjects ; do
    ## for subject in 335_A ; do	
	
	for timepoint in A C ; do
	    
	    subjectNumber=${subject%%_*}
	    dicomDir=$dicomRoot/${subjectNumber}_${timepoint}
	    if [[ -d ${dicomDir} ]] ; then
		sdir=$( cd ${dicomDir} ; ls -d s* | head -1 ) 
		iFile=$( cd ${dicomDir} ; find $sdir -regex  ".*/i.*\\(MRDC\\|CFMRI\\).*" -print -quit )
		dicomFile=${dicomDir}/${iFile}
		
		if [[ -f ${dicomFile}  ]] ; then

		    ## acquisition date
		    acqDate=$( dicom_hdr ${dicomFile} | grep "ID Acquisition Date" | awk -F "//" '{print $3}' )

		    echo "$subjectNumber,$group,$timepoint,$acqDate"
		    
		else
		    ## echo "*** Cannot read ${dicomFile}. Skipping."
		    echo "$subjectNumber,$group,$timepoint,NF"
		fi
		
		
	    else
		## echo "** *Cannot find ${dicomDir}. Skipping"
		echo "$subjectNumber,$group,$timepoint,ND"
		
	    fi

	done
	## exit	
    done
done
