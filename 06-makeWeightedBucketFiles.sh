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
GETOPT_OPTIONS=$( $GETOPT  -o "l:e:c" --longoptions "seedlist:,excessiveMotionThresholdPercentage:,cleaned" -n ${programName} -- "$@" )
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
	-e|--excessiveMotionThresholdPercentage)
	    excessiveMotionThresholdPercentage=$2; shift 2 ;;	
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

if [[ "x$excessiveMotionThresholdPercentage" == "x" ]] ; then
    echo "No excessiveMotionThresholdPercentage threshold was set. Exiting"
    exit 255
else
    echo "****************************************************************************************************"
    echo "*** Using ${excessiveMotionThresholdPercentage}% as motion cutoff threshold"
    echo "****************************************************************************************************"

fi

function makeFwhmFiles {
    local groupName="$1"
    local subjectList="$2"
    local groupData="$3"
    
    remlFwhmFile=$groupData/restingstate.$groupName.fwhmEstimates.tab
    cat /dev/null > $remlFwhmFile
    for subject in $subjectList ; do
	if [[ ! -f $DATA/$subject/rsfcPreprocessed/00_DO_NOT_ANALYSE_${subject}_${excessiveMotionThresholdPercentage}percent.txt ]] ; then
	    
	    cat $DATA/$subject/rsfcPreprocessed/${subject}.cleanEPI.blur.est.1D  >> $remlFwhmFile
	fi
    done
}


function makeAutocorrelatedBuckets {
    local groupName=$1
    local subjectList="$2"

    ## Autocorrelation based seeds
    for seed in $seeds ; do

	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi
	
	csvFilePrefix=subjectOrder.$groupName.${seedName}
	noDataCsvFile=$GROUP_DATA/nodata.$groupName.${seedName}.txt
	doNotAnalyzeCsvFile=$GROUP_DATA/doNotAnalyse.$groupName.${seedName}.txt

	cat /dev/null > $noDataCsvFile
	cat /dev/null > $doNotAnalyzeCsvFile

	echo "subject" > $GROUP_DATA/$csvFilePrefix.csv
	bucketListFile=$GROUP_DATA/bucketList.$groupName.${seedName}.txt
	rm -f $bucketListFile
	touch $bucketListFile

	problemSubjectList=""
	subbrikLabels=""
	for subject in $subjectList; do
	    zScoreFile=$DATA/$subject/rsfc/$seedName${cleanedSuffix}/${seedName}.z-score+tlrc.HEAD
	    if [[ ! -f $DATA/$subject/rsfcPreprocessed/00_DO_NOT_ANALYSE_${subject}_${excessiveMotionThresholdPercentage}percent.txt ]] ; then
		if [[ -f $zScoreFile ]] ; then
		    echo "$zScoreFile" >> $bucketListFile
		    echo "$subject" >> $GROUP_DATA/$csvFilePrefix.csv
		    if [[ -z $subbrikLabels ]] ; then
			subbrikLabels="$subject"			
		    else
			subbrikLabels="$subbrikLabels $subject"
		    fi
		else 
		    echo "$subject" >> $noDataCsvFile
		    echo "*** WARNING $zScoreFile does not exist!"
		fi
	    else 
		echo "$subject" >> $doNotAnalyzeCsvFile
		echo "*** WARNING: Found $DATA/$subject/rsfcPreprocessed/00_DO_NOT_ANALYSE_${subject}_${excessiveMotionThresholdPercentage}percent.txt. Skipping ${subject}."
	    fi
	done ## end of for subject in $subjects; do
	restingStateBucketFile=$GROUP_DATA/restingstate.bucket.$groupName.${seedName}+tlrc
	echo "Making resting state bucket file $restingStateBucketFile"

	rm -f ${restingStateBucketFile}.*
	3dbucket -fbuc -prefix    $restingStateBucketFile filelist:$bucketListFile
	3drefit  -relabel_all_str "$subbrikLabels" ${restingStateBucketFile}
	## rm -f $bucketListFile

	# rm -f ${restingStateBucketFile%%+*}.masked*
	# if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then
	#     echo "*** Using $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD as mask"
	#     3dcalc -datum float -a $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -b ${restingStateBucketFile} -expr "a*b" -prefix ${restingStateBucketFile%%+*}.masked
	# else 
	#     3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${restingStateBucketFile} -expr "a*b" -prefix ${restingStateBucketFile%%+*}.masked
	# fi
	# rm -f ${restingStateBucketFile}.*

    done ## end of for seed in $seeds ; do  

    if [[ ! -z $problemSubjectList ]] ; then 
	echo "*** The following subjects has non existant z-score files:"
	echo $problemSubjectList
    fi
}

function maskBuckets {
    local groupName=$1
    for seed in $seeds ; do

	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi
	
	restingStateBucketFile=$GROUP_DATA/restingstate.bucket.$groupName.${seedName}+tlrc
	rm -f ${restingStateBucketFile%%+*}.masked*
	if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then
	    echo "*** Using $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD as mask"
	    3dcalc -datum float -a $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -b ${restingStateBucketFile} -expr "a*b" -prefix ${restingStateBucketFile%%+*}.masked
	else 
	    3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${restingStateBucketFile} -expr "a*b" -prefix ${restingStateBucketFile%%+*}.masked
	fi
	rm -f ${restingStateBucketFile}.*
    done ## end of for seed in $seeds ; do 

}

function makeScaledDataLink {
    local dd="$1"
    ## rootname = rname
    local rname=$(dirname $dd)
    ## dirname = dname 
    local dname=$(basename $dd)
    
    local scaleddname=$( echo $dname | sed 's/diff/scaled.diff/' )

    if [[ -d $rname/$scaleddname ]] ; then 
	( cd $rname ; ln -sf $dname $scaleddname )
    fi
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
    GROUP_RESULTS=$DATA/Group.results${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/clean.ncl.subjectList.txt )"
    mddSubjects="$( cat ../data/config/clean.mdd.subjectList.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    echo "Calling makeGreyMatterMask.sh"
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddAndCtrl $seedName

    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
    
fi

## Group.data.BDI.diff.withAandC  Group.data.CDRSR.diff.withAandC  Group.data.MASC.tscore.diff.withAandC        Group.results.BDI.diff.withAandC  Group.results.CDRSR.diff.withAandC  Group.results.MASC.tscore.diff.withAandC
## Group.data.CDI.diff.withAandC  Group.data.CGAS.diff.withAandC   Group.data.RADS.Total.Tscore.diff.withAandC  Group.results.CDI.diff.withAandC  Group.results.CGAS.diff.withAandC   Group.results.RADS.Total.Tscore.diff.withAandC

if [[ $makeCdrsrbuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CDRS.t.score.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.CDRS.t.score.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.CDRS.t.score.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.CDRS.t.score.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"

    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA

    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA

    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly

fi

if [[ $makeBdiBuckets == 1 ]] ; then
    GROUP_DATA=$DATA/Group.data.BDI.II.Total.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.BDI.II.Total.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.BDI.II.Total.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.BDI.II.Total.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA
    
    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
fi

if [[ $makeCgasBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CGAS.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.CGAS.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.CGAS.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.CGAS.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA
    
    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
fi

if [[ $makeMascBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.MASC.tscore.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.MASC.tscore.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.MASC.tscore.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.MASC.tscore.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA
    
    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
fi

if [[ $makeCdiBuckets == 1 ]] ; then 
    GROUP_DATA=$DATA/Group.data.CDI.Total.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.CDI.Total.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.CDI.Total.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.CDI.Total.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA
    
    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
fi

if [[ $makeRadsBuckets == 1 ]] ; then  
    GROUP_DATA=$DATA/Group.data.RADS.Total.tscore.diff.withAandC${cleanedSuffix}
    GROUP_RESULTS=$DATA/Group.results.RADS.Total.tscore.diff.withAandC${cleanedSuffix}
    ctrlSubjects="$( cat ../data/config/new.ncl.subjectList.with.RADS.Total.tscore.AandC.txt )"
    mddSubjects="$( cat ../data/config/new.mdd.subjectList.with.RADS.Total.tscore.AandC.txt )"
    subjects="$ctrlSubjects $mddSubjects"
    
    [[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
    [[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS

    makeScaledDataLink $GROUP_DATA
    
    makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
    makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
    makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"

    makeFwhmFiles "mddAndCtrl" "$subjects"     $GROUP_DATA
    makeFwhmFiles "ctrlOnly"   "$ctrlSubjects" $GROUP_DATA
    makeFwhmFiles "mddOnly"    "$mddSubjects"  $GROUP_DATA
    
    ./makeGreyMatterMask.sh $GROUP_DATA $GROUP_RESULTS mddOnly $seedName
    maskBuckets mddAndCtrl
    maskBuckets ctrlOnly
    maskBuckets mddOnly
fi

#GROUP_DATA=$DATA/Group.data.withAandC
#GROUP_RESULTS=$DATA/Group.results.withAandC
#ctrlSubjects="$( cat ../data/config/ncl.subjectList.AandC.txt )"
#mddSubjects="$( cat ../data/config/mdd.subjectList.AandC.txt )"
#subjects="$ctrlSubjects $mddSubjects"
#[[ ! -d $GROUP_DATA ]]    && mkdir -p $GROUP_DATA
#[[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS
#makeAutocorrelatedBuckets "mddAndCtrl" "$subjects"
#makeAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
#makeAutocorrelatedBuckets "mddOnly"  "$mddSubjects"
