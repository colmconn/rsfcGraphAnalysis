#!/bin/bash

## set -x 

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



function fixDataTable {
    local dataTable="$1"

    lineCount=$( wc -l $dataTable | awk '{print $1}' )
    head -n $( expr $lineCount - 1 ) < $dataTable > ${dataTable}.new
    tail -1 $dataTable | sed "s/\\\\$//" >> ${dataTable}.new
    
    mv -f ${dataTable}.new ${dataTable}
}

cd $GROUP_RESULTS
grouping="mddAndCtrl"

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    
    dataTableFilename=${GROUP_DATA}/dataTable.${grouping}.${seedName}.txt
    fixDataTable $dataTableFilename
    if [[ -f $GROUP_RESULTS/mask.grey.$grouping.union.masked+tlrc.HEAD ]] ; then
	mask=$GROUP_RESULTS/mask.grey.$grouping.union.masked+tlrc.HEAD
    else
	mask=$MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz
    fi

    timestamp=$( date +%Y%m%d-%H%M%Z ) 

    echo "3dLME -prefix restingstate.mddAndCtrl.$seedName.lme.bucket.$timestamp \
	  -jobs 8 \
	  -mask $mask \
	  -qVars age.in.years \
	  -model'Group*Gender' \
	  -ranEff '~1' \
	  -SS_type 3 \
          -num_glt  2                                                   \
          -gltLabel 1 'MDD-NCL' -gltCode  1 'Group : 1*MDD -1*NCL'      \
          -gltLabel 2 'M-F'     -gltCode  2 'Gender : 1*M  -1*F'      \
	  -dataTable \
	  @$dataTableFilename"

done
    
