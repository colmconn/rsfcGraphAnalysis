#!/bin/bash

# set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log
GETOPT_OPTIONS=$( $GETOPT  -o "l:c" --longoptions "seedlist:,cleaned" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

cleaned=0
cleanedSuffix=""

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	-c|--cleaned)
	    cleaned=1; shift ;;	
	--) 
	    shift ; break ;;
	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ ! -f $seedList ] || [ "x$seedList" == "x" ] ; then
    echo "*** ERROR: The seed list file does not exit or was not provided. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList | sed "/#/d" ) )
fi

if [[ $cleaned -eq 1 ]] ; then
    echo "*** Will use cleaned versions of the data files"
    cleanedSuffix=".cleaned"

    GROUP_DATA=${GROUP_DATA}${cleanedSuffix}
    GROUP_RESULTS=${GROUP_RESULTS}${cleanedSuffix}
fi

cd $GROUP_RESULTS
groups="mddAndCtrl"

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi

    echo "*** Performing 1-sample t-test for the $seedName seed against zero (0)"
    if [[ -f mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then
	echo "*** Using mask.grey.$groups.union.masked+tlrc.HEAD as mask"
	3dttest++ -mask mask.grey.$groups.union.masked+tlrc.HEAD \
		  -setA $GROUP_DATA/restingstate.bucket.mddAndCtrl.${seedName}.masked+tlrc.HEAD  \
		  -prefix ttest.${seedName}
    else 
	3dttest++ -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz \
		  -setA $GROUP_DATA/restingstate.bucket.mddAndCtrl.${seedName}.masked+tlrc.HEAD  \
		  -prefix ttest.${seedName}
    fi
    
done
    
