#!/bin/bash

set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
#GROUP_DATA=$DATA/Group.data
#GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log
GETOPT_OPTIONS=$( $GETOPT  -o "l:is:n:op:c:" --longoptions "seedlist:,useInherentSmoothness,svc:,nn:,overwrite,pvalue:,cpvalue:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

useInherentSmoothness=0
overwrite=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	-p|--pvalue)
	    pValue=$2;
	    shift 2 ;;
	-c|--cpvalue)
	    cPvalue=$2;
	    shift 2 ;;
	-n|--nn)
	    NN=$2; 
	    shift 2 ;;
	-i|--useInherentSmoothness ) 
	    useInherentSmoothness=1; 
	    shift ;;
	-o|--overwrite ) 
	    overwrite=1; 
	    shift ;;
	--) 
	    shift ; break ;;
	
	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done



function makeBucketFilePrefix {
    local group=$1
    local seedName=$2
    local rvName=$3

    #       restingstate.mddOnly.R_SFA.3mm.MASC.tscore.rlm.bucket.20141015-1132PDT+tlrc.HEAD
    prefix="restingstate.${group}.${seedName}.${rvName}.rlm.bucket"

    echo $prefix
}

function pickLatestBucketFile {
    
    local prefix=$1
    local latest=$( ls -1t ${prefix}*+tlrc.HEAD | head -1 ) 

    if [ "x$latest" == "x" ] || [ ! -f $latest ] ; then
	exit
    fi
    echo $latest
}

function extractCoefBrikId {
    local rvName=$1
    local bucketFilename=$2
    
    label=$( 3dinfo -label $bucketFilename | tr "|" "\n" | grep "${rvName}.Value"  2> /dev/null )
    id=$( 3dinfo -label2index $label $bucketFilename 2> /dev/null )
    
    echo $id
}

function extractTvalueBrikId {
    local rvName=$1
    local bucketFilename=$2
    
    label=$( 3dinfo -label $bucketFilename | tr "|" "\n" | grep "${rvName}.t.value"  2> /dev/null )
    id=$( 3dinfo -label2index $label $bucketFilename 2> /dev/null )
    
    echo $id
}

function extractTStatpar {
    local bucket=$1
    local subbrikId=$2

    a=$(3dAttribute BRICK_STATSYM $bucket\[$subbrikId\] )
    b=${a##*(}
	c=${b%%)*}

    echo $c
}

function extractNVoxels {
    local nn=$1
    # corrected p-value
    local cPValue=$2
    ## voxelwise p-value
    local vPValue=$3
    local clustsimPrefix=$4

    #0.100  0.090  0.080  0.070  0.060  0.050  0.040  0.030  0.020  0.010
    # if [ $cPValue == 0.100 ] ; then
    # 	pvalueColumn=2
    # elif [ $cPValue == 0.090 ] ; then
    # 	pvalueColumn=3
    # elif [ $cPValue == 0.080 ] ; then
    # 	pvalueColumn=4
    # elif [ $cPValue == 0.070 ] ; then
    # 	pvalueColumn=5
    # elif [ $cPValue == 0.060 ] ; then
    # 	pvalueColumn=6
    # elif [ $cPValue == 0.050 ] ; then
    # 	pvalueColumn=7
    # elif [ $cPValue == 0.040 ] ; then
    # 	pvalueColumn=8
    # elif [ $cPValue == 0.030 ] ; then
    # 	pvalueColumn=9
    # elif [ $cPValue == 0.020 ] ; then
    # 	pvalueColumn=10
    # elif [ $cPValue == 0.010 ] ; then
    # 	pvalueColumn=11
    # fi


    if [[ ! -f $clustsimPrefix.NN$nn$side.1D ]] ; then
	echo "$clustsimPrefix.NN$nn.1D doesn't exist. Cannot continue, exiting."
	exit
    fi
    pvalueColumn=2
    nVoxels=$( cat $clustsimPrefix.NN$nn$side.1D | sed '/^#/d' | grep "^ $vPValue" | awk "{print \$$pvalueColumn}" )

    echo $nVoxels
}

if [ ! -f $seedList ] || [ "x$seedList" == "x" ] ; then
    echo "*** ERROR: The seed list file does not exit or was not provided. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList | sed "/#/d" ) )
fi

if [ "x$NN" == "x" ] ; then 
    ## nearest neighbour 1=touching at faces, 2=faces and edges 3=faces,
    ## edges and corners, just like in the afni clusterize window
    NN=1
fi
## the corrected p-value for the clusters. NB: the trailing zero at the
## end of this p-value is really important, do not omit it
##cPvalue=0.050
case $NN in
    1)
	rmm=1.01
	;;
    2)
	rmm=1.44
	;;
    3)
	rmm=1.75
	;;

    *) 
	echo "Unknown value ($NN) for NN. Exiting."
	exit 2 ;;
esac

if [ $useInherentSmoothness -eq 1 ] ; then 
    if [ -f $GROUP_DATA/${task}.mddAndCtrl.REML.fwhmEstimates.tab ] ; then 
	usedFwhm=$( cat $GROUP_DATA/${task}.mddAndCtrl.REML.fwhmEstimates.tab | gawk '{print $4}' | gawk '{s+=$1}END{print s/NR}' )
    else 
	usedFwhm=4.2
    fi
else 
    usedFwhm=4.2
fi

echo "*** Correcting for $usedFwhm mm smoothness"

if [[ "x$pValue" == "x" ]] ; then
    ## voxelwise pvalue
    pValue=0.05
    echo "*** Set voxelwise pvalue to $pValue (default)"
else
    echo "*** Set voxelwise pvalue to $pValue"
fi

if [[ "x$cPvalue" == "x" ]] ; then
    # clusterwise pvalue
    cPvalue=0.050
    echo "*** Set whole brain pvalue to $cPvalue (default)"	    
else
    useFirstColumn=1
    echo "*** Set whole brain pvalue to $cPvalue"    
fi

#[[ ! -d $GROUP_RESULTS ]] && mkdir -p $GROUP_RESULTS
#cd $GROUP_RESULTS

#regressionVariables="MASC.tscore"
#regressionVariables="MASC.tscore CDRS.tscore" # BDI.II"

regressionVariables="CDRSR.diff" # MASC.tscore.diff CGAS.diff RADS.Total.Tscore.diff"

## regressionVariables="BDI.diff"
#regressionVariables="CGAS.diff"
#regressionVariables="CDI.diff RADS.Total.Tscore.diff"

task="restingstate"

groups="mddOnly"

csvFile=regressions.parameters.csv

## for use with 3dClustSim
side="_1sided"

for regressionVariable in $regressionVariables ; do

    GROUP_DATA=$DATA/Group.data.$regressionVariable.withAandC
    GROUP_RESULTS=$DATA/Group.results.$regressionVariable.withAandC

    if [[ ! -d $GROUP_RESULTS ]] ; then
	echo "No such directory: $GROUP_RESULTS. Exiting"
	exit
    fi
    cd $GROUP_RESULTS

    if [[ $overwrite -eq 1 ]] || [[ ! -f $csvFile ]] ; then 
	echo "seed,polarity,usedFwhm,regressionVariable,coefBrikId,tvalueBrikId,rmm,nVoxels,df,tValue,pValue,cPvalue,nClusters,latestRlmBucketFile" > $csvFile
    fi

    echo "####################################################################################################"
    
    for seed in $seeds ; do
	
	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi
	
	## cstempPrefix=CStemp.fwhm${usedFwhm}
	cstempPrefix=CStemp.fwhm${usedFwhm}.pvalue.$pValue.cPvalue.$cPvalue
	if [ ! -f ${cstempPrefix}.NN1$side.1D ] ; then
	    echo "*** Running 3dClustSim"
	    #export OMP_NUM_THREADS=10
	    if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then 
		3dClustSim -mask mask.grey.$groups.union.masked+tlrc.HEAD -fwhm ${usedFwhm} -niml -both -prefix ${cstempPrefix} -pthr $pValue -athr $cPvalue
	    else 
		3dClustSim -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz  -fwhm ${usedFwhm} -niml -both -prefix ${cstempPrefix} -pthr $pValue -athr $cPvalue
	    fi
	fi
	
	bucketFilename=$GROUP_DATA/restingstate.bucket.$groups.${seedName}.masked+tlrc

	latestRlmBucketPrefix=$( makeBucketFilePrefix $groups $seedName $regressionVariable )
	latestRlmBucketFile=$( pickLatestBucketFile  $latestRlmBucketPrefix ) 
	echo "Group: $groups Regression Variable: $regressionVariable RLM bucket file: $latestRlmBucketFile"
	
	outputPrefix=${latestRlmBucketPrefix%%.rlm*}
	
	## extract positive and negative correlations
	
	coefBrikId=$( extractCoefBrikId $regressionVariable $latestRlmBucketFile )
	tvalueBrikId=$( extractTvalueBrikId $regressionVariable $latestRlmBucketFile )
	nVoxels=$( extractNVoxels $NN $cPvalue $pValue $cstempPrefix ) 
	df=$( extractTStatpar $latestRlmBucketFile $tvalueBrikId )
	tThreshold=$( cdf -p2t fitt $pValue $df | sed 's/t = //' )

	echo "### fwhm = ${usedFwhm}"
	echo "### coefBrikId = $coefBrikId"
	echo "### tvalueBrikId = $tvalueBrikId"
	echo "### rmm = $rmm"
	echo "### nVoxels = $nVoxels"
	echo "### df = $df"
	echo "### t-value = $tThreshold"
	echo "### voxelwise pValue = $pValue"
	echo "### corrected  pValue = $cPvalue"

	3dcalc -a $latestRlmBucketFile\[$coefBrikId\] -b $latestRlmBucketFile\[$tvalueBrikId\] -expr 'b*isnegative(a)' -prefix ${outputPrefix}.negative
	3dcalc -a $latestRlmBucketFile\[$coefBrikId\] -b $latestRlmBucketFile\[$tvalueBrikId\] -expr 'b*ispositive(a)' -prefix ${outputPrefix}.positive
	
	for polarity in positive negative ; do

	    suffix=regression.fwhm${usedFwhm}.$task.$groups.$seedName.and.$regressionVariable.$polarity
	    
	    # 3dclust -1Dformat -savemask clorder.regression.fwhm${usedFwhm}.$task.$groups.$regressionVariable.$polarity \
	    # 	-nosum -2thresh -$tThreshold $tThreshold -dxyz=1 $rmm $nVoxels -overwrite \
	    # 	${outputPrefix}.${polarity}+tlrc.HEAD  > clust.regression.fwhm${usedFwhm}.$task.$groups.$regressionVariable.$polarity.txt
	    
	    3dmerge -session . -prefix clorder.$suffix \
		    -1thresh $tThreshold \
		    -dxyz=1 -1clust_order $rmm $nVoxels \
		    -1erode 50 -1dilate \
		    ${outputPrefix}.${polarity}+tlrc.HEAD
	    3dclust -1Dformat -nosum -dxyz=1 $rmm $nVoxels clorder.$suffix+tlrc > clust.$suffix.txt
	    
	    # 3dclust -1Dformat -savemask clorder.$suffix \
	    # 	    -nosum -1thresh $tThreshold -dxyz=1 $rmm $nVoxels -overwrite \
	    # 	    ${outputPrefix}.${polarity}+tlrc.HEAD  > clust.$suffix.txt
	    
	    if [ -f clorder.$suffix+tlrc.HEAD ] ; then 
		3drefit -cmap INT_CMAP  clorder.$suffix+tlrc.HEAD

		3dcalc -a clorder.${suffix}+tlrc.HEAD -b ${outputPrefix}.${polarity}+tlrc.HEAD -expr "step(a)*b" -prefix clust.$suffix

		nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | awk '{print $1}' )
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${bucketFilename} > \
			   roiStats.$suffix.txt

		3dcalc -a $latestRlmBucketFile\[$coefBrikId\] -expr "a*is${polarity}(a)" -prefix rm.${outputPrefix}.coefficient.$polarity 
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD rm.${outputPrefix}.coefficient.$polarity+tlrc.HEAD   > roiStats.$suffix.averageCoefficientValue.txt
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${outputPrefix}.${polarity}+tlrc.HEAD                > roiStats.$suffix.averageTValue.txt
		rm -f rm.${outputPrefix}.coefficient.$polarity*
		
	    else
		nClusters=0
		echo "*** WARNING No clusters found!"
	    fi
	    echo "$seedName,$polarity,${usedFwhm},$regressionVariable,$coefBrikId,$tvalueBrikId,$rmm,$nVoxels,$df,$tThreshold,$pValue,$cPvalue,$nClusters,$latestRlmBucketFile" >> $csvFile
	    
	done ## end of for polarity in "positive negative" ; do

	##exit

    done ## end of for seed in $seeds ; do 

    cd $scriptsDir

    if [[ $overwrite -eq 0 ]] ; then 
	##echo "*** Making cluster location tables using Maximum intensity"
	##./cluster2Table.pl --space=mni --force -mi $GROUP_RESULTS

	echo "*** Making cluster location tables using Center of Mass"
	./cluster2Table.pl --space=mni $GROUP_RESULTS
    elif [[ $overwrite -eq 1 ]] ; then 
	##echo "*** Making cluster location tables using Maximum intensity"
	##./cluster2Table.pl --space=mni --force -mi $GROUP_RESULTS

	echo "*** Making cluster location tables using Center of Mass"
	./cluster2Table.pl --space=mni --force $GROUP_RESULTS
    fi
done ## end of for regressionVariable in $regressionVariables ; do



