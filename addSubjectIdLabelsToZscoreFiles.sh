#!/bin/bash

## set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data

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

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else
    #subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.nat.txt )"
    
    ## change task above too to match this list of subjects
    subjects=$( cd $DATA/; \ls -1d *_[ABC] )
fi

echo "*** Processing subjects: $subjects"

if [[ ! -f $seedList ]] || [[ "x$seedList" == "x" ]] ; then
    echo "*** ERROR: The seed list file does not exit or was not provided. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList | sed "/#/d" ) )
fi

if [[ $cleaned -eq 1 ]] ; then
    echo "*** Will use cleaned versions of the z-score files for each seed"
    cleanedSuffix=".cleaned"
fi

function addSubjectIdLabel {
    local subject=$1
    local zscoreFile="$2"
  
    if [[ -f $zscoreFile ]] ; then
	3drefit -sublabel 0 $subject $zscoreFile
    else
	echo "*** No such file: $zscoreFile"
    fi
}

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    
    for subject in $subjects ; do
	# zscoreFile=$DATA/${subject}/rsfc/${seedName}/${seedName}.z-score+tlrc.HEAD
	# addSubjectIdLabel $subject $zscoreFile

	zscoreFile=$DATA/${subject}/rsfc/${seedName}/${seedName}.z-score.masked+tlrc.HEAD
	addSubjectIdLabel $subject $zscoreFile

	
    done
done
