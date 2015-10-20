#!/bin/bash

#set -x 

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
GETOPT_OPTIONS=$( $GETOPT  -o "l:is:n:" --longoptions "seedlist:,useInherentSmoothness,svc:,nn:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

useInherentSmoothness=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	-n|--nn)
	    NN=$2; 
	    shift 2 ;;
	-i|--useInherentSmoothness ) 
	    useInherentSmoothness=1; 
	    shift ;;
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

if [ $useInherentSmoothness -eq 1 ] ; then 
    if [ -f $GROUP_DATA/${task}.mddAndCtrl.REML.fwhmEstimates.tab ] ; then 
	usedFwhm=$( cat $GROUP_DATA/${task}.mddAndCtrl.REML.fwhmEstimates.tab | gawk '{print $4}' | gawk '{s+=$1}END{print s/NR}' )
    else 
	usedFwhm=4.2
    fi
else 
    usedFwhm=4.2
fi

setA="MASC.tscore"

#setB="CDRS.tscore BDI.II"
#setB="BDI.II"

setB="RADS.Total.Tscore"

task="restingstate"

groups="mddOnly"

function makeClorderFilename {

    local seedName=$1
    local regressionVariable=$2
    local polarity=$3

    echo "clorder.regression.fwhm${usedFwhm}.$task.$groups.$seedName.and.$regressionVariable.$polarity+tlrc.HEAD"
    
}

cd $GROUP_RESULTS

for seed in $seeds ; do
    
    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    
    for set1 in $setA ; do
	
	for set1Polarity in positive negative ; do

	    set1Filename=$( makeClorderFilename $seedName $set1 $set1Polarity )

	    if [[ -f $set1Filename ]] ;then 
		for set2 in $setB ; do
		    for set2Polarity in positive negative ; do		
		    
			set2Filename=$( makeClorderFilename $seedName $set2 $set2Polarity )
			if [[ -f $set2Filename ]] ; then 
			    echo "####################################################################################################"
			    echo "set1: $set1Filename set2: $set2Filename"


			    ## compute the overlap

			    intersectionPrefix=intersection.$seedName.$set1.$set1Polarity.and.$set2.$set2Polarity
			    3dcalc -a $set1Filename -b $set2Filename -prefix $intersectionPrefix -expr "and(a, b)"
			    
			    ## compute what's unique to $set1Filename

			    set1OnlyPrefix=set1Only.$seedName.$set1.$set1Polarity.and.$set2.$set2Polarity
			    3dcalc -a $set1Filename -b $set2Filename -prefix temp -expr "step(a) + 2 * step(b)"
			    3dcalc -a temp+tlrc -prefix $set1OnlyPrefix -expr "equals(a, 1)"
			    rm -f temp+tlrc*

			    ## compute what's unique to $set2Filename
			    set2OnlyPrefix=set2Only.$seedName.$set1.$set1Polarity.and.$set2.$set2Polarity
			    3dcalc -a $set1Filename -b $set2Filename -prefix temp -expr "step(a) + 2 * step(b)"
			    3dcalc -a temp+tlrc -prefix $set2OnlyPrefix -expr "equals(a, 2)"
			    rm -f temp+tlrc*

			    ## now combine them so that we can see them in one image
			    3dcalc -a ${set1OnlyPrefix}+tlrc -b ${intersectionPrefix}+tlrc -c ${set2OnlyPrefix}+tlrc -expr "a + 3*b + 2*c" -prefix combined.$seedName.$set1.$set1Polarity.and.$set2.$set2Polarity
			    
			fi
		    done
		done
	    fi
	done
    done
done

