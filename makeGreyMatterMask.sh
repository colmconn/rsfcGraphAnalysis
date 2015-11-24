#!/bin/bash

## set -x

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard

GETOPT_OPTIONS=$( $GETOPT  -o "d:r:s:g:j:" --longoptions "data:,results:,seed:,groups:,subjects:" -n ${programName} -- "$@" )
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
	-d|--data)
	    group_data=$2; shift 2
	    echo "*** group_data is in $group_data"
	    ;;
	-r|--results)
	    group_results=$2; shift 2
	    echo "*** group_results is in $group_results"
	    ;;
	-s|--seed)
	    seedName=$2; shift 2
	    echo "*** seedName is $seedName"
	    ;;
	-g|--groups)
	    groups=$2; shift 2
	    echo "*** groups is $groups"
	    ;;
	-j|--subjects)
	    subjectsFile=$2
	    if [[ -f $subjectsFile ]] ; then
		nsubjects=$( cat $subjectsFile | sed 1d | wc -l ) 
		echo "*** Reading $nsubjects subjects from $subjectsFile" 
		subjects="$subjects $( cat $subjectsFile | sed 1d )"
	    else
		echo "*** No such file: $subjectsFile"
	    fi
	    shift 2
	    ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $group_data ] ; then 
    echo "*** ERROR: The group data directory was not provided."
    exit
fi
if [ -z $group_results ] ; then 
    echo "*** ERROR: The group results directory was not provided."
    exit
fi
if [ -z $seedName ] ; then 
    echo "*** ERROR: The seed name was not provided."
    exit
fi
if [ -z $groups ] ; then 
    echo "*** ERROR: The groups name was not provided."
    exit
fi
if [ -z "$subjects" ] ; then 
    echo "*** No subjects list was provided. Will use the list in the subjectOrder file"
fi


## exit

##regressionVariables="CDRSR.diff MASC.tscore.diff CGAS.diff RADS.Total.Tscore.diff"

## regressionVariables="BDI.diff"
## regressionVariables="CGAS.diff"
## regressionVariables="CDI.diff RADS.Total.Tscore.diff"

task="restingstate"

## groups="mddOnly"


function makeMask {

    if [[ ! -d $group_data ]] ; then
	echo "No such directory: $group_data. Exiting"
	exit
    fi

    if [[ ! -d $group_results ]] ; then
	echo "No such directory: $group_results. Exiting"
	exit
    fi

    if [[ -z "$subjects" ]] ; then
	echo "*** Loading subjects list from $group_data/$subjectOrderFilename"
	## seedName=L_BLA.3mm
	subjectOrderFilename=subjectOrder.$groups.${seedName}.csv
	
	subjects=$( cat $group_data/$subjectOrderFilename | sed 1d )
    fi
    
   
    gmList=""
    for subject in $subjects ; do
	
	cd $DATA/$subject/rsfcPreprocessed/tmp
	
	if [[ ! -f ../$subject.mask.grey.MNI.nii.gz ]] ; then 

	    echo "*** $subject"
	    ## now warp the GM to MNI space for this subject
	    3dbucket -prefix mask.grey.nii mask.grey+orig.HEAD

	    3dresample -master ../../anat/${subject}.anat_struc.std.nii.gz -prefix ./mask.grey.resampled.nii \
    		       -inset mask.grey+orig.HEAD

	    flirt -in mask.grey.nii.gz -ref mask.grey.resampled.nii.gz \
	      -out mask.grey.func2anat.flirt -omat mask.grey.func2anat.flirt.mat

	    applywarp \
		--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
		--in=mask.grey.nii.gz \
		--warp=../../anat/${subject}.std.2.MNI.warpcoef \
		--premat=mask.grey.func2anat.flirt.mat \
		--interp=nn \
		--out=mask.grey.MNI
	    
	    mv -f mask.grey.MNI.nii.gz ../$subject.mask.grey.MNI.nii.gz
	fi
	gmList="$gmList $subject/rsfcPreprocessed/$subject.mask.grey.MNI.nii.gz"
	
    done
    echo "The following grey matter masks are going into group masks:"
    echo $gmList
    
    cd $DATA
    3dMean -mask_inter -prefix mask.grey.$groups.intersection $gmList
    3dMean -mask_union -prefix mask.grey.$groups.union $gmList

}


# if [[ $# -gt 0 ]] ; then
#     GROUP_DATA="$1"    ## $DATA/Group.data
#     GROUP_RESULTS="$2" ## $DATA/Group.results
#     groups="$3"        ## "mddOnly"
#     seedName="$4"      ## "L_BLA.3mm"
    
## makeMask $GROUP_DATA $GROUP_RESULTS $groups $seedName

makeMask
## pwd
## ls -l mask*
mv -f $DATA/mask.grey.$groups.* $group_results

for ff in $group_results/mask.grey.$groups.*HEAD ; do
    if [[ ! -f ${ff%%+*}.masked+tlrc.HEAD ]] ; then
	echo "*** Masking with MNI152 T1 3mm brain mask"
	3dcalc -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b $ff -expr "a*b" -prefix ${ff%%+*}.masked
    fi
done

# else

#     echo "You provided too few comand line arguments. Cannot continue."
#     exit
    
#     # GROUP_DATA=$DATA/Group.data
#     # GROUP_RESULTS=$DATA/Group.results
#     # groups="mddAndCtrl"
#     # seedName="$4"
    
#     # makeMask $GROUP_DATA $GROUP_RESULTS $groups
#     # mv -f $DATA/mask.grey.$groups.* $GROUP_RESULTS

#     # for ff in $GROUP_RESULTS/mask.grey.$groups.*HEAD ; do
#     # 	3dcalc -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b $ff -expr "a*b" -prefix ${ff%%+*}.masked 
#     # done

#     # groups="mddOnly"
#     # for regressionVariable in $regressionVariables ; do
	
#     # 	GROUP_DATA=$DATA/Group.data.$regressionVariable.withAandC
#     # 	GROUP_RESULTS=$DATA/Group.results.$regressionVariable.withAandC

#     # 	makeMask $GROUP_DATA $GROUP_RESULTS $groups
#     # 	mv -f $DATA/mask.grey.$groups.* $GROUP_RESULTS

#     # 	for ff in $GROUP_RESULTS/mask.grey.$groups.*HEAD ; do
#     # 	    3dcalc -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b $ff -expr "a*b" -prefix ${ff%%+*}.masked 
#     # 	done
#     # done
# fi
