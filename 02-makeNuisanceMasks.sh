#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

GETOPT_OPTIONS=$( $GETOPT  -o "s:f" --longoptions "subject:,force" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

## force recomputation of masks? 1=yes 0=no
force=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-f|--force)
	    echo "*** Forcing recomputation of masks"
	    force=1; shift ;;
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

cd $DATA/$subjectNumber 

if [[ -f anat/segment_seg_0.nii.gz ]] && [[ -f anat/segment_seg_1.nii.gz ]] && [[ -f anat/segment_seg_2.nii.gz ]] ; then 

    cd anat
    
    ln -sf segment_seg_0.nii.gz csf.nii.gz
    ln -sf segment_seg_1.nii.gz gm.nii.gz
    ln -sf segment_seg_2.nii.gz wm.nii.gz

    if [[ $force -eq 1 ]] ; then
	echo "*** Deleting outputs to force recomputation of transformation matrix and masks"
	rm -f \
	    MNI.affine2anat_struc_brain.std.mat HO-sub-maxprob-thr50-native.nii* \
	    HO-ventricles.native.nii* HO-wm.native.nii* ventricles.mask* white.mask* gm+orig.*
    fi

    ## invert the affine transform from native T1 space to MNI to yield a MNI -> native T1 transform
    [[ ! -f MNI.affine2anat_struc_brain.std.mat ]] && convert_xfm \
	-omat MNI.affine2anat_struc_brain.std.mat \
	-inverse ${subjectNumber}.anat_struc_brain.std.2.MNI.affine.mat 

    ## use the transform created above to warp the HO subcortical atlas into native space
    [[ ! -f HO-sub-maxprob-thr50-native.nii.gz ]] &&  flirt \
	-ref ${subjectNumber}.anat_struc_brain.std  \
	-in $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz \
	-out HO-sub-maxprob-thr50-native \
	-applyxfm \
	-init MNI.affine2anat_struc_brain.std.mat \
	-interp nearestneighbour

    ## now extract the ventricles and white matter
    [[ ! -f HO-ventricles.native.nii.gz ]] && 3dcalc -a HO-sub-maxprob-thr50-native.nii.gz -expr "or(equals(a, 3), equals(a, 14))" -prefix HO-ventricles.native.nii
    [[ ! -f HO-wm.native.nii.gz ]]         && 3dcalc -a HO-sub-maxprob-thr50-native.nii.gz -expr "or(equals(a, 1), equals(a, 12))" -prefix HO-wm.native.nii

    [[ ! -f ventricles.mask+orig.HEAD ]] && 3dcalc \
	-a csf.nii.gz \
	-b HO-ventricles.native.nii.gz \
	-c ${subjectNumber}.anat_struc_brain.mask.std.nii.gz \
	-expr "step(a) * step(b) * step(c)" -prefix ventricles.mask
    
    [[ ! -f white.mask+orig.HEAD ]] && 3dcalc \
	-a wm.nii.gz  \
	-b HO-wm.native.nii.gz         \
	-c ${subjectNumber}.anat_struc_brain.mask.std.nii.gz \
	-expr "step(a) * step(b) * step(c)" -prefix white.mask

    [[ ! -f gm+orig.HEAD ]] && 3dcopy gm.nii.gz gm

else
    echo "*** FAST must have been run to create the $(pwd)/anat/segment_seg_[123].nii.gz files before this script can be run."
    echo "*** Run 01-preprocessAnatomy.sh first."
    echo "*** Exiting."

fi
