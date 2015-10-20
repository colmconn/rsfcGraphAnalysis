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
GETOPT_OPTIONS=$( $GETOPT  -o "s:l:" --longoptions "subject:,seedlist:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided. Exiting"
    exit
fi

if [ ! -f $seedList ] ; then
    echo "*** ERROR: The seed list file does not exit. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList ) )
fi

echo "*** Computing RSFC for the following seeds:"
echo $seeds

[[ ! -d $DATA/$subjectNumber/rsfc ]] && mkdir $DATA/$subjectNumber/rsfc
cd $DATA/$subjectNumber/rsfc

preprocessedRsfcDir=$DATA/$subjectNumber/rsfcPreprocessed

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi

    mkdir ${seedName}

    echo "*** Extracting timeseries for seed ${seed}"
    3dROIstats -quiet -mask_f2short -mask ${seed} ${preprocessedRsfcDir}/${subjectNumber}.pm.cleanEPI.MNI.nii.gz > ${seedName}/${seedName}.ts.1D

    echo "*** Computing Correlation for seed ${seedName}"
    3dfim+ -input ${preprocessedRsfcDir}/${subjectNumber}.pm.cleanEPI.MNI.nii.gz -ideal_file ${seedName}/${seedName}.ts.1D -out Correlation -bucket ${seedName}/${seedName}_corr
    
    echo "*** Z-transforming correlations for seed ${seedName}"
    3dcalc -datum float -a ${seedName}/${seedName}_corr+tlrc.HEAD -expr 'log((a+1)/(a-1))/2' -prefix ${seedName}/${seedName}.z-score

    3drefit -sublabel 0 $subjectNumber ${seedName}/${seedName}.z-score+tlrc.HEAD
done
