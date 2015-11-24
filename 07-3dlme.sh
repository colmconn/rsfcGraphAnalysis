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

if [ ! -f $seedList ] || [ "x$seedList" == "x" ] ; then
    echo "*** ERROR: The seed list file does not exit or was not provided. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList | sed "/#/d" ) )
fi

if [[ $cleaned -eq 1 ]] ; then
    echo "*** Will use cleaned versions of the data files"
    cleanedSuffix=".cleaned"

    GROUP_DATA=${GROUP_DATA}${cleanedSuffix}
    GROUP_RESULTS=${GROUP_RESULTS}${cleanedSuffix}
fi



function fixDataTable {
    local dataTable="$1"
    echo "Fixing $dataTable"

    lineCount=$( wc -l $dataTable | awk '{print $1}' )
    head -n $( expr $lineCount - 1 ) < $dataTable > ${dataTable}.new
    tail -1 $dataTable | sed 's/\\$//' >> ${dataTable}.new

    mv -f ${dataTable}.new ${dataTable}

}

cd $GROUP_RESULTS
grouping="mddAndCtrl"

taskFile="$scriptsDir/3dlme-taskFile.txt"
cat /dev/null > $taskFile

for seed in $seeds ; do

    seedName=${seed##*/}
    if echo $seedName | grep -q "nii" ; then 
	seedName=${seedName%%.nii*}
    else 
	seedName=${seedName%%+*}
    fi
    
    if [[ -f $GROUP_RESULTS/mask.grey.$grouping.union.masked+tlrc.HEAD ]] ; then
	mask=$GROUP_RESULTS/mask.grey.$grouping.union.masked+tlrc.HEAD
    else
	mask=$MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz
    fi


    dataTableFilename=${GROUP_DATA}/dataTable.${grouping}.group.and.gender.${seedName}.txt
    fixDataTable $dataTableFilename
    timestamp=$( date +%Y%m%d-%H%M%Z ) 

    ## ###############################################################
    ## NOTE WELL NOTE WELL NOTE WELL NOTE WELL NOTE WELL NOTE WELL
    ## ###############################################################
    ## When creating contrasts to investigate interations like in GLTs
    ## 3 and 4 below, ensure that the constant level (e.g., MDD in GLT
    ## 3 and NCL in GLT 4) is always on the left. If you do not, the
    ## graphLmeRois.r script will not correctly graph your results for
    ## you.
    ## ###############################################################
    ## NOTE WELL NOTE WELL NOTE WELL NOTE WELL NOTE WELL NOTE WELL
    ## ###############################################################

    if [[ 1 == 1 ]] ; then 
    cmd="3dLME      -prefix $GROUP_RESULTS/restingstate.mddAndCtrl.group.and.gender.$seedName.lme.bucket.$timestamp \
	  -jobs     8 \
	  -mask     $mask \
	  -model    'Group*Gender+Full' \
	  -ranEff   '~1' \
          -qVars    'Full' \
	  -SS_type  3 \
          -num_glt  4 \
          -gltLabel 1 'MDD-NCL'     -gltCode 1 'Group  : 1*MDD -1*NCL' \
          -gltLabel 2 'M-F'         -gltCode 2 'Gender : 1*M   -1*F' \
          -gltLabel 3 'M.MDD-M.NCL' -gltCode 3 'Gender : 1*M Group : 1*MDD -1*NCL' \
          -gltLabel 4 'F.MDD-F.NCL' -gltCode 4 'Gender : 1*F Group : 1*MDD -1*NCL' \
          -num_glf  2 \
          -glfLabel 1 'MDD-NCL.GLF' -glfCode 1 'Gender : 1*M & 1*F     Group  : 1*MDD -1*NCL' \
          -glfLabel 2 'M-F.GLF'     -glfCode 2 'Group  : 1*MDD & 1*NCL Gender : 1*M -1*F' \
	  -dataTable \
	  @$dataTableFilename"
    #echo "3dLME command is: $cmd"
    echo "$cmd" >> $taskFile
    fi
    
    dataTableFilename=${GROUP_DATA}/dataTable.${grouping}.group.and.age.${seedName}.txt
    fixDataTable $dataTableFilename
    # timestamp=$( date +%Y%m%d-%H%M%Z ) 

    if [[ 0 == 1 ]] ; then 
    cmd="3dLME      -prefix $GROUP_RESULTS/restingstate.mddAndCtrl.group.and.age.$seedName.lme.bucket.$timestamp \
    	  -jobs     8 \
    	  -mask     $mask \
    	  -model    'Group*age+Full' \
    	  -ranEff   '~1' \
          -qVars    'age,Full' \
    	  -SS_type  3 \
          -num_glt  4                                                 \
          -gltLabel 1 'MDD-NCL'    -gltCode  1 'Group : 1*MDD -1*NCL' \
          -gltLabel 2 'age'        -gltCode  2 'age :' \
          -gltLabel 3 'MDD-AgeEff' -gltCode  3 'Group : 1*MDD age :'  \
          -gltLabel 4 'NCL-AgeEff' -gltCode  4 'Group : 1*NCL age :'  \
    	  -dataTable \
    	  @$dataTableFilename"
    #echo "3dLME command is: $cmd"
    echo "$cmd" >> $taskFile
    fi

done
    
## jobname
#$ -N rsfcLme

## queue
#$ -q all.q

## binary?
#$ -b y

## rerunnable?
#$ -r y

## merge stdout and stderr?
#$ -j y

## send no mail
#$ -m n

## execute from the current working directory
#$ -cwd

## use a shell to run the command
#$ -shell yes 
## set the shell
#$ -S /bin/bash

## preserve environment
#$ -V 

nTasks=$( cat $taskFile | wc -l )

## run the enqueuing command (qsub) in a subshell to avoid replacing the
## current process witht he qsub command which is what happend when
## using the exec command. This way we can continue to run commands
## (e.g., qstat) in the parent shell after the qsub command has been
## executed
(
sge_command="qsub -N amygdalaRsfcLme -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
#echo $sge_command
echo "Queuing job... "
exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF
)

qstat
