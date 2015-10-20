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

logDir=${DATA}/log
GETOPT_OPTIONS=$( $GETOPT  -o "l:" --longoptions "seedlist:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ ! -f $seedList ] ; then
    echo "*** ERROR: The seed list file does not exit. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList ) )
fi

echo "*** Extracting GM information for the following seeds:"
echo $seeds

if [[ ! -d $DATA/vbm/stats ]] ; then
    echo "No such directory $DATA/vbm/stats. Fun the VBM scripts first."
    exit
fi

cd $DATA/vbm/stats
[[ ! -d rois ]] && mkdir rois

if [[ -f GM_mod_merg_s2.nii.gz ]] ; then 
    for seed in $seeds ; do

	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi

	echo "*** Extracting timeseries for seed ${seed}"
	3dROIstats -quiet -mask_f2short -mask ${seed} GM_mod_merg_s2.nii.gz > rois/${seedName}.gm.1D

    done
fi
