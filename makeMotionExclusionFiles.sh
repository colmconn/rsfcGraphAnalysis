#!/bin/bash

## set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard

subjects=$( cd $DATA ; ls -d [0-9][0-9][0-9]_[ABC]* )
#subjects=$*

motionThresholdPrecentage=0.2
threshold=20

for subject in $subjects ; do

    if [[ ! -d $DATA/$subject/rsfcPreprocessed ]] ; then
	continue
    fi
    cd $DATA/$subject/rsfcPreprocessed
    
    if [[ -f ${subject}.pm.cleanEPI+orig.HEAD ]] ; then 

	numberOfCensoredVolumes=$( cat ${subject}.pm.censor.1D | gawk '{a+=(1-$0)}END{print a}' )
	totalNumberOfVolumes=$( cat ${subject}.pm.censor.1D | wc -l )
	cutoff=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

	if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then
	    rm -f 00_DO_NOT_ANALYSE*
	    echo "*** A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than $threshold % (n=$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_${subject}_${threshold}percent.txt
	    echo "*** WARNING: $subject will not be analysed due to having more that $threshold % of their volumes censored."
	fi
    fi
done
