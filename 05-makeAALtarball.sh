#!/bin/bash

set -x

trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data

GETOPT_OPTIONS=$( $GETOPT  -o "l:e:" --longoptions "excessiveMotionThresholdPercentage:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-e|--excessiveMotionThresholdPercentage)
	    excessiveMotionThresholdPercentage=$2; shift 2 ;;	
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [[ "x$excessiveMotionThresholdPercentage" == "x" ]] ; then
    echo "No excessiveMotionThresholdPercentage threshold was set. Exiting"
    exit 255
else
    echo "****************************************************************************************************"
    echo "*** Using ${excessiveMotionThresholdPercentage}% as motion cutoff threshold"
    echo "****************************************************************************************************"

fi

function makeDirectoryList {
    subjects="$1"
    task="$2"
    directoryListFile="$3"

    # echo "$subjects"
    # echo "$task"
    # echo "$directoryListFile"


    (( i=0 ))
    (( e=0 ))
    (( n=0 ))        
    for subject in $subjects ; do
cleanedEpiFile=$DATA/$subject/${task}Preprocessed/${subject}.pm.cleanEPI.MNI.nii.gz
	if [[ ! -f $DATA/$subject/${task}Preprocessed/00_DO_NOT_ANALYSE_${subject}_${excessiveMotionThresholdPercentage}percent.txt ]] ; then
	    if [[ -d $DATA/$subject/${task}NoBP/yeo7lib ]] ; then
		echo $subject/${task}NoBP/yeo7lib/ >> $directoryListFile
#	echo "$i, Include, $subject"
		(( i=i+1 ))
	    else
		echo "$n, ****No data*******, $subject"
		(( n=n+1 ))		
	    fi
	else
 echo "$e, *******Exclude*******, $subject"
	    (( e=e+1 ))			    
	fi
    done
    echo "Total included: $i"
    echo "Total excluded: $e"
    echo "Total nodata  : $n"
}

#subjects=$( cd $DATA; \ls -1d *_[ABC] )
#subjects="364_A 364_C 370_A 421_A"
#subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/estop.subjList.txt )"
## subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/rs114a.subj.txt )"
subjects="$( cat ../data/config/flex.rsfc.mddAndCtrl.csv )"
#task=rsfc
task=estop
cat /dev/null > rm.${task}.yeo7splitEstopNoBP.txt
makeDirectoryList "$subjects" ${task} rm.${task}.yeo7splitEstopNoBP.txt
echo -n "Making tarball..."
(cd $DATA/; tar czf ${task}.yeo7split.tar.gz $( cat ../scripts/rm.${task}.yeo7splitEstopNoBP.txt ) )
echo "Done"

# subjects="$( cat /data/sanDiego/cESTOP/data/config/clean.estop.subjList.txt /data/sanDiego/cESTOP/data/config/clean.estop.mddfollowup.txt  )"
# task=estop 
# cat /dev/null > rm.${task}.aalList.txt
# makeDirectoryList "$subjects" ${task} rm.${task}.aalList.txt
# (cd $DATA; tar czf ${task}.aal.tar.gz $( cat ../scripts/rm.${task}.aalList.txt ) )
