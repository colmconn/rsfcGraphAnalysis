#!/bin/bash

## set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
##GROUP_DATA=$DATA/Group.data
##GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log
GETOPT_OPTIONS=$( $GETOPT  -o "l:is:n:op:c:s:evd:r:" \
			   --longoptions "seedlist:,useInherentSmoothness,svc:,nn:,overwrite,pvalue:,cpvalue:,sided,cleaned,covaried,data:,results:" \
			   -n ${programName} -- "$@" )
exitStatus=$?
if [[ $exitStatus != 0 ]] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

useInherentSmoothness=0
overwrite=0

cleaned=0
cleanedSuffix=""

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
	-s|--sided )
	    ss=$2
	    if [[ $ss == "1sided" ]] ; then 
		side="1sided"
	    elif [[ $ss == "2sided" ]] ; then 
		side="2sided"
	    elif [[ $ss == "bisided" ]] ; then 
		side="bisided"
	    else
		echo "Unknown argument provided to -s or --sided. Valid values are 1sided, 2sided, bisided. Defaulting to 1sided"
		side="_1sided"		
	    fi
	    shift 2 ;;	
	-o|--overwrite ) 
	    overwrite=1; 
	    shift ;;
	-e|--cleaned)
	    cleaned=1; shift ;;		
	-v|--covaried)
	    covaried=1; shift ;;
	-d|--data)
	    GROUP_DATA=$2;
	    shift 2 ;;
	-r|--results)
	    GROUP_RESULTS=$2;
	    shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

function extractTStatpars {
    local bucket="$1"
    local subbrikId="$2"

    a=$(3dAttribute BRICK_STATSYM $bucket"[$subbrikId]" )
    b=${a##*(}
    c=${b%%)*}

    echo $( echo $c | tr "," " " )
}

function makeClustSimFilePrefix {
    local group=$1
    local fwhm="$2"

    local pvalue="$3"
    local cpvalue="$4"

    if [[ -z "$pvalue" ]] || [[ -z "$cpvalue" ]] ; then 
	prefix="CStemp.fwhm${usedFwhm}.${group}"
    else
	prefix=CStemp.fwhm${usedFwhm}.pvalue.$pvalue.cPvalue.$cpvalue
    fi
    echo $prefix
}

function extractNVoxels {
    local nn=$1
    # corrected p-value
    local cPValue=$2
    ## voxelwise p-value
    local vPValue=$3
    local clustsimPrefix=$4
    local side="$5"

    #0.100  0.090  0.080  0.070  0.060  0.050  0.040  0.030  0.020  0.010
    # if using the default whole brain p values used by 3dClustSim
    # then the following if/elsif/fi statement should be uncommented

    pValueColumn=$( cat $clustsimPrefix.NN${nn}_${side}.1D | grep "#[ ]*pthr" | sed "s/#\s*pthr\s*|\s*//" | tr ' ' '\n' | grep -nh "\\.05"  | awk -F ':' '{print $1}' )


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

    ## pvalueColumn=2
    if [[ "X$pValueColumn" == "X" ]] ; then
	nVoxels="NA"
    else
	pValueColumn=$( expr $pValueColumn + 1 )
	nVoxels=$( cat $clustsimPrefix.NN${nn}_${side}.1D | sed '/^#/d' | grep "^ $vPValue" | awk "{print \$$pValueColumn}" )
    fi
    
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

if [[ $cleaned -eq 1 ]] ; then
    echo "*** Will use cleaned versions of the data files"
    cleanedSuffix=".cleaned"

    GROUP_DATA=${GROUP_DATA}${cleanedSuffix}
    GROUP_RESULTS=${GROUP_RESULTS}${cleanedSuffix}
fi

if [[ $covaried -eq 1 ]] ; then
    echo "*** Will use covaried versions of the data files"
    covariedSuffix=".covaried"
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

if [[ "x$side" == "x" ]] ; then
    echo "*** No value provided for side. Defaulting to 1sided"
    side="_1sided"
else
    echo "*** Running a $side test"
fi

if [[ "x$GROUP_DATA" == "x" ]] ; then
    echo "*** No value provided for GROUP_DATA (-d or --data). Cannot continue."
    exit
fi

if [[ "x$GROUP_RESULTS" == "x" ]] ; then
    echo "*** No value provided for GROUP_RESULTS (-r or --results). Cannot continue."
    exit
fi

GROUP_DATA=$( readlink -f $GROUP_DATA )
if [[ ! -d "$GROUP_DATA" ]] ; then
    echo "*** No such directory: $GROUP_DATA"
    echo "Cannot continue."
    exit
fi

GROUP_RESULTS=$( readlink -f $GROUP_RESULTS )
if [[ ! -d "$GROUP_RESULTS" ]] ; then
    echo "*** No such directory: $GROUP_RESULTS"
    echo "Cannot continue."
    exit
fi

echo "*** Will use data          files in $GROUP_DATA"
echo "*** Will use group results files in $GROUP_RESULTS"
cd $GROUP_RESULTS

groups="mddAndCtrl"
csvFile=parameters.fwhm${usedFwhm}.$groups.csv

if [[ $overwrite -eq 1 ]] || [[ ! -f $csvFile ]] ; then 
    echo "analysis,seed,fwhm,tLabel,tValueBrikId,tThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,bucketFile" > $csvFile
fi

tLabelPrefix="MDD-NCL"

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi

    ## cstempPrefix=CStemp.fwhm${usedFwhm}.pvalue.$pValue.cPvalue.$cPvalue
    ## cstempPrefix=$( makeClustSimFilePrefix $groups $usedFwhm $pValue $cPvalue )
    cstempPrefix="CStemp.fwhmxyz7.98x8.05x7.32" #.NN1_1sided.1D 
    # if [ ! -f ${cstempPrefix}.NN${NN}_${side}.1D ] ; then
    # 	echo "*** Running 3dClustSim"
    # 	echo "*** Output will be saved to files begining with: $cstempPrefix"
    # 	#export OMP_NUM_THREADS=10
    # 	if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then 
    # 	    3dClustSim -mask mask.grey.$groups.union.masked+tlrc.HEAD -fwhm ${usedFwhm} -both -prefix ${cstempPrefix} -pthr $pValue -athr $cPvalue
    # 	else 
    # 	    3dClustSim -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz  -fwhm ${usedFwhm} -both -prefix ${cstempPrefix} -pthr $pValue -athr $cPvalue
    # 	fi
    # 	mv -f 3dClustSim.cmd ${cstempPrefix}.3dClustSim.cmd
    # fi

    bucketFilename=$GROUP_DATA/restingstate.bucket.$groups.${seedName}.masked+tlrc
    tTestFile=ttest.${groups}.${seedName}${covariedSuffix}+tlrc
    
    #addStatTableCmd="$( cat ${cstempPrefix}.3dClustSim.cmd ) $tTestFile"
    #echo "*** Statistic table addition command is: $addStatTableCmd"
    #eval "$addStatTableCmd"
    
    tValueBrikId=$( 3dinfo -label2index "${tLabelPrefix}_Tstat" $tTestFile 2> /dev/null )
    tContrastBrikId=$( 3dinfo -label2index "${tLabelPrefix}_mean" $tTestFile 2> /dev/null )    
    ## tContrastBrikId=$( expr $tValueBrikId - 1 )

    nVoxels=$( extractNVoxels $NN $cPvalue $pValue ${cstempPrefix} $side )
    if [[ "$nVoxels" == "NA" ]] ; then
	echo "*** Couldn't get the correct number of voxels to go with pvalue=$pValue and corrected pvalue=$cPvalue"
	echo "*** You may need to pad these values with zeros to ensure you match the correct row and column in $cstempPrefix.NN${NN}_${side}.1D"
	exit
    fi
    df=$( extractTStatpars "$tTestFile" "${tLabelPrefix}_Tstat" )

    tThreshold=$( cdf -p2t fitt $pValue $df | sed 's/t = //' )
    echo "### fwhm = ${usedFwhm}"
    echo "### tLabelPrefix = $tLabelPrefix"
    echo "### tContrastBrikId = $tContrastBrikId"
    echo "### tValueBrikId = $tValueBrikId"
    echo "### tThreshold = $tThreshold"
    echo "### rmm = $rmm"
    echo "### nVoxels = $nVoxels"
    echo "### df = $df"
    echo "### voxelwise pValue = $pValue"
    echo "### corrected  pValue = $cPvalue"

    suffix=fwhm${usedFwhm}.$groups.$seedName
    mddOnlySuffix=fwhm${usedFwhm}.mddOnly.$seedName
    ctrlOnlySuffix=fwhm${usedFwhm}.ctrlOnly.$seedName

    3dmerge -session . -prefix clorder.$suffix \
	    -2thresh -$tThreshold $tThreshold \
	    -1clust_order $rmm $nVoxels \
	    -dxyz=1 \
	    -1erode 50 -1dilate \
	    -1dindex $tContrastBrikId -1tindex $tValueBrikId  \
	    $tTestFile
    
    3dclust -1Dformat -nosum -dxyz=1 $rmm $nVoxels clorder.$suffix+tlrc.HEAD > clust.$suffix.txt
    nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )
    ## if [[ -f clorder.$suffix+tlrc.HEAD ]] ; then
    if [[ $nClusters -gt 0 ]] ; then 

	3dcalc -a clorder.${suffix}+tlrc.HEAD -b ${tTestFile}\[$tValueBrikId\] -expr "step(a)*b" -prefix clust.$suffix

	3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename         > roiStats.$suffix.txt
	#3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $ctrlOnlyBucketFilename > roiStats.$ctrlOnlySuffix.txt
	#3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $mddOnlyBucketFilename  > roiStats.$mddOnlySuffix.txt

	3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD         ${tTestFile}\[$tContrastBrikId\] > roiStats.$suffix.averageContrastValue.txt
	3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD         ${tTestFile}\[$tValueBrikId\]    > roiStats.$suffix.averageTValue.txt
	3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD -minmax ${tTestFile}\[$tValueBrikId\]    > roiStats.$suffix.minmaxTValue.txt
	
	echo "$df" > text.$suffix.degreesOfFreedom.txt
	3drefit -cmap INT_CMAP clorder.$suffix+tlrc.HEAD
    else
	## nClusters=0
	echo "*** WARNING No clusters found!"
    fi
    echo "NA,$seedName,${usedFwhm},$tLabel,$tValueBrikId,$tThreshold,$rmm,$nVoxels,$df,$pValue,$cPvalue,$nClusters,$tTestFile" >> $csvFile
done

cd $scriptsDir
##echo "*** Making cluster location tables using Maximum intensity"
##./cluster2Table.pl --space=mni --force -mi $GROUP_RESULTS

#echo "*** Making cluster location tables using Center of Mass"
#./cluster2Table.pl --space=mni --force $GROUP_RESULTS
