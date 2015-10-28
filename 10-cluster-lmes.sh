#!/bin/bash

## set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts
logDir=${DATA}/log

task="restingstate"

groups="mddAndCtrl"

function makeBucketFilePrefix {
    local group=$1
    local seedName=$2
    local analysis="$3"

    ## restingstate",        groups, seedName, "lme.bucket",       sep=".")
    if [[ ! -z "$analysis" ]] ;then 
	prefix="restingstate.${group}.${analysis}.${seedName}.lme.bucket"
    else 
	prefix="restingstate.${group}.${seedName}.lme.bucket"
    fi
    echo $prefix
}

function makeClustSimFilePrefix {
    local group=$1
    local fwhm="$2"

    prefix="CStemp.fwhm${usedFwhm}.${group}"

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

function extractFStatpars {
    local bucket=$1
    local subbrikId=$2

    a=$(3dAttribute BRICK_STATSYM $bucket\[$subbrikId\] )
    b=${a##*(}
    c=${b%%)*}

    echo $( echo $c | tr "," " " )
}

function extractTStatpars {
    local bucket=$1
    local subbrikId=$2

    a=$(3dAttribute BRICK_STATSYM $bucket\[$subbrikId\] )
    b=${a##*(}
    
    echo ${b%%)*}
}


function fixStatLabel {

    echo "$( echo $1 | awk 'BEGIN  { OFS="."}; {print $1,$2}' | sed 's/:/X/g' )"
    # echo "$( echo $1 | sed "s/:/X/g" )"
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
    if [ $cPValue == 0.100 ] ; then
	pvalueColumn=2
    elif [ $cPValue == 0.090 ] ; then
	pvalueColumn=3
    elif [ $cPValue == 0.080 ] ; then
	pvalueColumn=4
    elif [ $cPValue == 0.070 ] ; then
	pvalueColumn=5
    elif [ $cPValue == 0.060 ] ; then
	pvalueColumn=6
    elif [ $cPValue == 0.050 ] ; then
	pvalueColumn=7
    elif [ $cPValue == 0.040 ] ; then
	pvalueColumn=8
    elif [ $cPValue == 0.030 ] ; then
	pvalueColumn=9
    elif [ $cPValue == 0.020 ] ; then
	pvalueColumn=10
    elif [ $cPValue == 0.010 ] ; then
	pvalueColumn=11
    fi

    if [[ "X$pvalueColumn" == "X" ]] ; then
	nVoxels="NA"
    else
	nVoxels=$( cat $clustsimPrefix.NN${nn}_${side}.1D | sed '/^#/d' | grep "^ $vPValue" | awk "{print \$$pvalueColumn}" )
    fi

    echo $nVoxels
}

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

fi

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

#export OMP_NUM_THREADS=8

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

cd $GROUP_RESULTS

analysis="group.and.gender" 
csvFile=parameters.fwhm${usedFwhm}.$groups.csv
echo "analysis,seed,fwhm,f/zLabel,f/zValueBrikId,f/zThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,bucketFile" > $csvFile

for seed in $seeds ; do
	
    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    
    echo "####################################################################################################"
    echo "### Seed is: $seedName"
    
    bucketFilename=$GROUP_DATA/$task.bucket.$groups.${seedName}.masked+tlrc.HEAD
    ctrlOnlyBucketFilename=$GROUP_DATA/$task.bucket.ctrlOnly.${seedName}.masked+tlrc.HEAD
    mddOnlyBucketFilename=$GROUP_DATA/$task.bucket.mddOnly.${seedName}.masked+tlrc.HEAD
    
    if [[ ! -f $bucketFilename ]] ; then
	pwd
	echo "*** ERROR cannot find $bucketFilename. Exiting"
	exit
    fi
    
    latestLmeBucketPrefix=$( makeBucketFilePrefix $groups $seedName  $analysis )
    latestLmeBucketFile=$( pickLatestBucketFile  $latestLmeBucketPrefix ) 
    echo "### Most recent bucket file is $latestLmeBucketFile"

    OIFS="$IFS"
    IFS=';'
    fValueBrikLabels=( $( 3dinfo -label $latestLmeBucketFile | tr "|" "\n" | grep " F$" | grep -v "Intercept" | tr "\n" ";" 2> /dev/null ) )
    zValueBrikLabels=( $( 3dinfo -label $latestLmeBucketFile | tr "|" "\n" | grep "Z$" | grep -v "Intercept" | tr "\n" ";" 2> /dev/null ) )
    IFS="$OIFS"
    
    echo "*** Got the following F value brik labels from the LME bucket: ${fValueBrikLabels[@]}"
    echo "*** Got the following Z value brik labels from the LME bucket: ${zValueBrikLabels[@]}"

    side="1sided"
    
    cstempPrefix=$( makeClustSimFilePrefix $groups $usedFwhm )
    if [[ ! -f ${cstempPrefix}.NN${NN}_${side}.1D ]] ; then
	echo "*** Running 3dClustSim"
	echo "*** Output will be saved to files begining with: $cstempPrefix"
	export OMP_NUM_THREADS=40
	if [[ -f $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD ]] ; then 
	    3dClustSim -mask $GROUP_RESULTS/mask.grey.$groups.union.masked+tlrc.HEAD -fwhm ${usedFwhm} -both -prefix ${cstempPrefix}
	else 
	    3dClustSim -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz  -fwhm ${usedFwhm} -both -prefix ${cstempPrefix}
	fi

	mv -f 3dClustSim.cmd ${cstempPrefix}.3dClustSim.cmd
    fi
    addStatTableCmd="$( cat ${cstempPrefix}.3dClustSim.cmd ) $latestLmeBucketFile"
    echo "*** Statistic table addition command is: $addStatTableCmd"
    eval "$addStatTableCmd"

    ## ####################################################################################################
    ## Clustering begins in ernest here
    ## ####################################################################################################    
    
    for fLabel in "${fValueBrikLabels[@]}" ; do
	
	fixedFLabel=$( fixStatLabel "$fLabel" ) 
	fValueBrikId=$( 3dinfo -label2index "$fLabel" $latestLmeBucketFile 2> /dev/null )
	nVoxels=$( extractNVoxels $NN $cPvalue $pValue ${cstempPrefix} $side )
	if [[ "$nVoxels" == "NA" ]] ; then
	    echo "*** Couldn't get the correct number of voxels to go with pvalue=$pValue and corrected pvalue=$cPvalue"
	    echo "*** You may need to pad these values with zeros to ensure you match the correct row and column in $cstempPrefix.NN${NN}_${side}.1D"
	    exit
	fi
	df=$( extractFStatpars $latestLmeBucketFile $fValueBrikId )
	fThreshold=$( cdf -p2t fift $pValue $df | sed 's/t = //' )
	
	echo "### **************************************"
	echo "### Processing the \"$fLabel\" subbrik"
	echo "### Parameters are:"
	echo "### $latestLmeBucketFile: \'$fLabel\' $fValueBrikId"
	#echo "### fwhm = ${fwhm[3]}"
	echo "### fwhm = ${usedFwhm}"
	echo "### fLabel = $fLabel"
	echo "### fixedFLabel = $fixedFLabel"
	echo "### fValueBrikId = $fValueBrikId"
	echo "### fThreshold = $fThreshold"
	echo "### rmm = $rmm"
	echo "### nVoxels = $nVoxels"
	echo "### df = $df"
	echo "### voxelwise pValue = $pValue"
	echo "### corrected  pValue = $cPvalue"

	suffix=fwhm${usedFwhm}.$task.$groups.$seedName.$fixedFLabel
	#mddOnlySuffix=fwhm${usedFwhm}.$task.mddOnly.$seedName.$fixedFLabel
	#ctrlOnlySuffix=fwhm${usedFwhm}.$task.ctrlOnly.$seedName.$fixedFLabel
	
	# #3dclust -1Dformat -savemask clorder.$suffix -nosum -1dindex $fValueBrikId -1tindex $fValueBrikId -1noneg -2thresh -$fThreshold $fThreshold \
	# #    -dxyz=1 $rmm $nVoxels $latestLmeBucketFile  > clust.$suffix.txt

	3dmerge -session . -prefix clorder.$suffix -1noneg -2thresh -$fThreshold $fThreshold -dxyz=1 -1clust_order $rmm $nVoxels -1erode 50 -1dilate -1dindex $fValueBrikId -1tindex $fValueBrikId $latestLmeBucketFile
	3dclust -1Dformat -nosum -dxyz=1 $rmm $nVoxels clorder.$suffix+tlrc.HEAD > clust.$suffix.txt

	if [ -f clorder.$suffix+tlrc.HEAD ] ; then 

	    3dcalc -a clorder.${suffix}+tlrc.HEAD -b ${latestLmeBucketFile}\[$fValueBrikId\] -expr "step(a)*b" -prefix clust.$suffix

	    nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )

	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename         > roiStats.$suffix.txt
	    #3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $ctrlOnlyBucketFilename > roiStats.$ctrlOnlySuffix.txt
	    #3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $mddOnlyBucketFilename  > roiStats.$mddOnlySuffix.txt
	    
	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${latestLmeBucketFile}\[$fValueBrikId\]         > roiStats.$suffix.averageFvalue.txt

	     3drefit -cmap INT_CMAP clorder.$suffix+tlrc.HEAD
	else
	    nClusters=0
	    echo "*** WARNING No clusters found!"
	fi
	echo "$analysis,$seedName,${usedFwhm},$fLabel,$fValueBrikId,$fThreshold,$rmm,$nVoxels,$df,$pValue,$cPvalue,$nClusters,$latestLmeBucketFile" >> $csvFile
    done ## end of for fLabel in $fValueBrikLabels ; do

    echo "####################################################################################################"

    side="2sided"
    for zLabel in "${zValueBrikLabels[@]}" ; do

	fixedZLabel=$( fixStatLabel "$zLabel" ) 
	zValueBrikId=$( 3dinfo -label2index "$zLabel" $latestLmeBucketFile 2> /dev/null )
	zContrastBrikId=$( expr $zValueBrikId - 1 )
	nVoxels=$( extractNVoxels $NN $cPvalue $pValue ${cstempPrefix} $side )
	if [[ "$nVoxels" == "NA" ]] ; then
	    echo "*** Couldn't get the correct number of voxels to go with pvalue=$pValue and corrected pvalue=$cPvalue"
	    echo "*** You may need to pad these values with zeros to ensure you match the correct row and column in $cstempPrefix.NN${NN}_${side}.1D"
	    exit
	fi
	## df=$( extractTStatpars $latestLmeBucketFile $tValueBrikId )
	zThreshold=$( cdf -p2t fizt $pValue | sed 's/t = //' )

	echo "### **************************************"
	echo "### Processing the \"$zLabel\" subbrik"
	echo "### Parameters are:"
	echo "### $latestLmeBucketFile: \'$zLabel\' $zValueBrikId"
	#echo "### fwhm = ${fwhm[3]}"
	echo "### fwhm = ${usedFwhm}"
	echo "### zLabel = $zLabel"
	echo "### fixedZLabel = $fixedZLabel"
	echo "### zContrastBrikId = $zContrastBrikId"
	echo "### zValueBrikId = $zValueBrikId"
	echo "### zThreshold = $zThreshold"
	echo "### rmm = $rmm"
	echo "### nVoxels = $nVoxels"
	## echo "### df = $df"
	echo "### voxelwise pValue = $pValue"
	echo "### corrected  pValue = $cPvalue"

	suffix=fwhm${usedFwhm}.$task.$groups.$seedName.$fixedZLabel
	# #mddOnlySuffix=fwhm${usedFwhm}.$task.mddOnly.$seed.$fixedTLabel
	# #ctrlOnlySuffix=fwhm${usedFwhm}.$task.ctrlOnly.$seed.$fixedTLabel

	3dmerge -session . -prefix clorder.$suffix -2thresh -$zThreshold $zThreshold -dxyz=1 -1clust_order $rmm $nVoxels -1erode 50 -1dilate -1dindex $tContrastBrikId -1tindex $tValueBrikId  $latestLmeBucketFile  
	3dclust -1Dformat -nosum -dxyz=1 $rmm $nVoxels clorder.$suffix+tlrc.HEAD > clust.$suffix.txt

	if [[ -f clorder.$suffix+tlrc.HEAD ]] ; then 

	    3dcalc -a clorder.${suffix}+tlrc.HEAD -b ${latestLmeBucketFile}\[$zContrastBrikId\] -expr "step(a)*b" -prefix clust.$suffix
  
	    nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )
	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename         > roiStats.$suffix.txt
	    #3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $ctrlOnlyBucketFilename > roiStats.$ctrlOnlySuffix.txt
	    #3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $mddOnlyBucketFilename  > roiStats.$mddOnlySuffix.txt

	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${latestLmeBucketFile}\[$zContrastBrikId\]      > roiStats.$suffix.averageContrastValue.txt
	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${latestLmeBucketFile}\[$zValueBrikId\]         > roiStats.$suffix.averageZValue.txt

	    3drefit -cmap INT_CMAP clorder.$suffix+tlrc.HEAD
	else
	    nClusters=0
	    echo "*** WARNING No clusters found!"
	fi
	echo "$analysis,$seedName,${usedFwhm},$zLabel,$zValueBrikId,$zThreshold,$rmm,$nVoxels,NA,$pValue,$cPvalue,$nClusters,$latestLmeBucketFile" >> $csvFile
    done

done ## end of for seed in $seeds ; do

cd $scriptsDir
##echo "*** Making cluster location tables using Maximum intensity"
##./cluster2Table.pl --space=mni --force -mi $GROUP_RESULTS

echo "*** Making cluster location tables using Center of Mass"
./cluster2Table.pl --space=mni $GROUP_RESULTS
