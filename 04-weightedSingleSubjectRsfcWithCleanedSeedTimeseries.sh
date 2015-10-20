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
    seeds=( $( eval echo $( cat $seedList ) ) ) 
fi

echo "*** Computing RSFC for the following seeds:"
echo ${seeds[*]}

[[ ! -d $DATA/$subjectNumber/rsfc ]] && mkdir $DATA/$subjectNumber/rsfc
rsfcDir=$DATA/$subjectNumber/rsfc
cd $rsfcDir

preprocessedRsfcDir=$DATA/$subjectNumber/rsfcPreprocessed

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed


function getSeedName {
    local seed=$1
    
    seedName=$( basename ${seed} ) 
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    echo $seedName
}

function getOtherSeedIndices {

    local unwantedSeedIndex=$1

    declare -a wantedSeedIndices

    ((xx=0))
    for (( yy=0; ${yy} < $numberOfSeeds; yy++ ))
    do
	if [[ $yy -ne $unwantedSeedIndex ]] ; then
	    wantedSeedIndices[${xx}]=$yy
	    ((xx++))
	fi
    done
    echo "${wantedSeedIndices[@]}"
}

# get length of an array
numberOfSeeds=${#seeds[@]}

if [[ $numberOfSeeds -lt 2 ]] ; then
    echo "*** The number of seeds must be at least 2 to use this script. You have $numberOfSeeds seeds in your seed file"
    exit
fi
 
# use for loop read all nameservers
for (( ii=0; ii<${numberOfSeeds}; ii++ ));
do
    seed=${seeds[$ii]}
    seedName=$( getSeedName ${seed} )

    seedDir=${seedName}.cleaned
    if [[ -d ${seedDir} ]] ; then
	rm -rf ${seedDir}
    fi
    mkdir ${seedDir}

    for (( jj=0; jj<${numberOfSeeds}; jj++ ));
    do
    	sd=${seeds[$jj]}
    	sn=$( getSeedName ${sd} )
    
    	echo "*** Extracting timeseries for seed ${sd}"
    	#echo "*** ii=$ii, seed=$seed, seedname=$seedName, jj=$jj, sd=$sd, sn=$sn"

	3dcalc -datum float -a ${sd} -b ${preprocessedRsfcDir}/${subjectNumber}.pm.cleanEPI.MNI.nii.gz -expr "a*b" -prefix ${preprocessedRsfcDir}/__tmp__${seedName}__${subjectNumber}.pm.cleanEPI.MNI
	3dmaskave -quiet -mask ${sd} ${preprocessedRsfcDir}/__tmp__${seedName}__${subjectNumber}.pm.cleanEPI.MNI+tlrc.HEAD  > ${seedDir}/${sn}.ts.1D
	rm -f ${preprocessedRsfcDir}/__tmp__${seedName}__${subjectNumber}.pm.cleanEPI.MNI+tlrc.*
	echo "***"
    	
    done
    # echo "*** ii=$ii, seed=$seed, seedname=$seedName"
    # getOtherSeedIndices $ii

    otherSeedIndices=( $( getOtherSeedIndices $ii ) ) 
    # echo "*** unwantedIndex=$ii, otherSeedIndices = ${otherSeedIndices[@]}"

    ## now try building the command line to clean the target timeseries
    cleaningCommand="$scriptsDir/cleanSeedTimeseries.r -G --session $rsfcDir/${seedDir} --prefix ${seedName}.cleaned --LHS $( getSeedName $seed ).ts.1D --RHS \""
    for (( zz=0 ; zz<${#otherSeedIndices[@]}; zz++ ))
    do
	cleaningCommand="${cleaningCommand}$( getSeedName ${seeds[${otherSeedIndices[$zz]}]} ).ts.1D "
    done
    cleaningCommand="$cleaningCommand\""
    
    echo "*** Cleaning command: $cleaningCommand"

    eval "$cleaningCommand"
    
    echo "*** Computing Correlation for seed ${seedName}"
    3dfim+ -input ${preprocessedRsfcDir}/${subjectNumber}.pm.cleanEPI.MNI.nii.gz -ideal_file ${seedDir}/${seedName}.cleaned.ts.1D -out Correlation -bucket ${seedDir}/${seedName}_corr
    
    echo "*** Z-transforming correlations for seed ${seedName}"
    3dcalc -datum float -a ${seedDir}/${seedName}_corr+tlrc.HEAD -expr 'log((a+1)/(a-1))/2' -prefix ${seedDir}/${seedName}.z-score 

done
