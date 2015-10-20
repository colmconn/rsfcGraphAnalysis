#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=/data/sanDiego/rsfcGraphAnalysis
DATA=$ROOT/data
FREESURFER_SUBJECTS_DIR=/data/sanDiego/freesurferAnalysis/data
scriptsDir=${ROOT}/scripts

function convertMgzToAfni {
    echo "*** Converting and gzipping $1 to $2"
    mri_convert -it mgz -ot nii -i $1 $2
    gzip -9 $2
}


GETOPT_OPTIONS=$( $GETOPT  -o "s:dk" --longoptions "subject:,destrieux,killiany" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

atlas="Destrieux"
aseg="aparc.a2009s+aseg.mgz"

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-d|--destrieux)
	    atlas="Destrieux"
	    aseg="aparc.a2009s+aseg.mgz"; shift 1 ;;
	-k|--killiany)
	    atlas="Desikan-Killiany"
	    aseg="aparc+aseg.mgz"; shift 1 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided."
    exit
fi

echo "*** Using the $atlas segmentation file: $aseg" 

cd $DATA/$subjectNumber
asegPrefix=${aseg%%.mgz}

if  [[ -f $FREESURFER_SUBJECTS_DIR/$subjectNumber/mri/$aseg ]] && \
    [[ -f $DATA/$subjectNumber/$subjectNumber.pm.cleanEPI+orig.HEAD  ]] && \
    [[ -f $DATA/$subjectNumber/$subjectNumber.pm.cleanEPI+orig.BRIK.gz  ]] ; then 

    [[ ! -f ${asegPrefix}.nii.gz ]] && convertMgzToAfni $FREESURFER_SUBJECTS_DIR/$subjectNumber/mri/$aseg ${asegPrefix}.nii

    if [[ -f ${asegPrefix}.nii.gz ]] ; then 

	[[ ! -f ${asegPrefix}.2.cleanEPI+orig.HEAD ]] && \
	    3dresample -rmode NN \
	    -master $subjectNumber.pm.cleanEPI+orig \
	    -prefix ${asegPrefix}.2.cleanEPI \
	    -inset  ${asegPrefix}.nii.gz
	
	echo "*** Extracting timeseries from each ROI in the mask"
	3dROIstats -quiet -mask_f2short \
	    -mask ${asegPrefix}.2.cleanEPI+orig \
	    $subjectNumber.pm.cleanEPI+orig > ${subjectNumber}.rsfc.1D
    else
	echo "*** Conversion of $aseg to NIfTI falied for subject $subjectNumber. Cannot continue."
	exit 1
    fi

    
fi
