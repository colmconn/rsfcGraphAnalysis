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

## GROUP_DATA=$DATA/Group.data.acc.amygdala.paper
## GROUP_RESULTS=$DATA/Group.results.acc.amygdala.paper


## GROUP_DATA=$DATA/Group.data.acc.graph.paper
## GROUP_RESULTS=$DATA/Group.results.acc.graph.paper


## GROUP_DATA=$DATA/Group.data.kaiser.amygdala.paper
## GROUP_RESULTS=$DATA/Group.results.kaiser.amygdala.paper


## GROUP_DATA=$DATA/Group.data.kaiser.graph.paper
## GROUP_RESULTS=$DATA/Group.results.kaiser.graph.paper





MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

# taskFile=taskFile-mask-z-score.txt
# cat /dev/null > $taskFile
# for bucketListFile in $GROUP_DATA/bucketList.mddAndCtrl.* ; do
#     groups="$( echo ${bucketListFile##*/} | cut -d '.' -f 2 )"
#     if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then 
# 	for zscoreFile in $( cat $bucketListFile ) ; do
# 	    echo "3dcalc -datum float -a $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -b ${zscoreFile} -expr \"a*b\" -prefix ${zscoreFile%%+*}.masked" >> $taskFile
# 	done
#     else
# 	echo "*** No such file $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD. Skipping."
#     fi
# done

# parallel --wd $( pwd ) -j16 < $taskFile


taskFile=taskFile-mask-z-score.txt
#subjects="320_A 323_A 149_A 150_A 161_A 300_A 309_A 316_A 338_A 344_A 345_A 349_A 361_A"
## seedList="../data/config/juelich_amygdala_seeds_weights.txt "

##seedList="../data/config/juelich_whole_amygdala_seeds.txt"

## seedList="../data/config/kaiser_supplemental_seeds.txt"
## seedList="../data/config/dlpfc_seeds.txt"

## seedList="../data/config/short_ACC_seed_list.txt"

# seedList="../data/config/juelich_whole_amygdala_seeds.txt
# ../data/config/short_ACC_seed_list.txt
# ../data/config/hippocampus_ventral_striatum_seeds.txt
# ../data/config/followup-dlpfc-ins-IP-MPFC-seeds.txt
# ../data/config/Fox-Goldapple-seeds.txt"

## seedList="../data/config/dmn_sn_yeo7network_seeds.txt"

## seedList="../data/config/miller-dmn.txt"
## seedList="../data/config/jacobs-seeds.txt"
## seedList="../data/config/goldapple-ofc-seeds.txt"
## seedList="../data/config/gabbay-striatum-seeds.txt"

## seedList="../data/config/tremblay-seeds.txt"
## seedList="../data/config/goldapple-vlpfc-seeds.txt"
seedList="../data/config/goldapple-dlpfc-seeds.txt"

seeds=$( eval echo $( cat ${seedList} | sed "/#/d" ) )
groups="mddAndCtrl"

echo "Masking RSFC z-score fores for the following seeds:"
echo $seeds
## exit

cat /dev/null > $taskFile

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi

    ## subjects=$( cat $GROUP_DATA/subjectOrder.mddAndCtrl.${seedName}.csv | sed 1d )
    ## subjects=$( cat $GROUP_DATA/subjectOrder.mddAndCtrl.L_whole_amygdala.3mm.csv | sed 1d )

    subjects=$( cd ../data/ ; ls -d [0-9][0-9][0-9]_{A,B,C,D,E} [0-9][0-9][0-9]_{A,B,C,D,E}2 2> /dev/null )
    for ss in $subjects ; do
	ddir=/data/sanDiego/rsfcGraphAnalysis/data/$ss/rsfc/$seedName
	zscoreFile=$ddir/${seedName}.z-score+tlrc.HEAD
	## ( cd $ddir ; 3dcalc -datum float -a $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -b ${zscoreFile} -expr "a*b" -prefix ${zscoreFile%%+*}.masked )
	echo "3dcalc -datum float -a $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -b ${zscoreFile} -expr \"a*b\" -prefix ${zscoreFile%%+*}.masked" >> $taskFile
    done
done
		 
parallel --wd $( pwd ) -j50 < $taskFile
