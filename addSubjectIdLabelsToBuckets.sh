#!/bin/bash

# set -x 

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
    local groupName=$1
    local groupData="$2"
    
    for seed in $seeds ; do

	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi
	
	csvFile=$groupData/subjectOrder.$groupName.${seedName}.csv
	subjectIds="$( cat $csvFile | sed 1d | tr "\n" " " )"

	echo "*** Adding the following subjects IDs as sub-brik labels"
	echo "*** $subjectIds"
	
	restingStateBucketFile=$groupData/restingstate.bucket.$groupName.${seedName}+tlrc
	maskedRestingStateBucketFile=${restingStateBucketFile%%+*}.masked+tlrc

	if [[ -f ${restingStateBucketFile}.HEAD ]] ; then 
	    3drefit  -relabel_all_str "$subjectIds" ${restingStateBucketFile}
	fi
	if [[ -f ${maskedRestingStateBucketFile}.HEAD ]] ; then 
	    3drefit  -relabel_all_str "$subjectIds" ${maskedRestingStateBucketFile}
	fi
    done
}

makeBetweenGroupBuckets=1
makeCdrsrbuckets=0
makeBdiBuckets=0
makeCgasBuckets=0
makeMascBuckets=0
makeCdiBuckets=0
makeRadsBuckets=0
seedName="L_BLA.weight.3mm"

if [[ $makeBetweenGroupBuckets == 1 ]] ; then 

    GROUP_DATA=$DATA/Group.data${cleanedSuffix}
    
    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

## Group.data.BDI.diff.withAandC  Group.data.CDRSR.diff.withAandC  Group.data.MASC.tscore.diff.withAandC        Group.results.BDI.diff.withAandC  Group.results.CDRSR.diff.withAandC  Group.results.MASC.tscore.diff.withAandC
## Group.data.CDI.diff.withAandC  Group.data.CGAS.diff.withAandC   Group.data.RADS.Total.Tscore.diff.withAandC  Group.results.CDI.diff.withAandC  Group.results.CGAS.diff.withAandC   Group.results.RADS.Total.Tscore.diff.withAandC

if [[ $makeCdrsrbuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CDRS.t.score.diff.withAandC${cleanedSuffix}

    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

if [[ $makeBdiBuckets == 1 ]] ; then
    GROUP_DATA=$DATA/Group.data.BDI.II.Total.diff.withAandC${cleanedSuffix}

    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

if [[ $makeCgasBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CGAS.diff.withAandC${cleanedSuffix}

    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

if [[ $makeMascBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.MASC.tscore.diff.withAandC${cleanedSuffix}

    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

if [[ $makeCdiBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CDI.Total.diff.withAandC${cleanedSuffix}

    if [[ ! -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

if [[ $makeRadsBuckets == 1 ]] ; then  
    GROUP_DATA=$DATA/Group.data.RADS.Total.tscore.diff.withAandC${cleanedSuffix}

    if [[ -d $GROUP_DATA ]] ; then
	addSubjectIdLabel "mddAndCtrl" $GROUP_DATA
	addSubjectIdLabel "ctrlOnly" $GROUP_DATA
	addSubjectIdLabel "mddOnly" $GROUP_DATA
    fi
fi

