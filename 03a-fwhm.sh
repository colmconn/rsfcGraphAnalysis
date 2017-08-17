#!/bin/bash

set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

subjects=$( cd $DATA ;  ls -d *_[ABC] )

taskFile=run/acf-est-TaskFile
cat /dev/null > $taskFile

for subject in $subjects ; do

    if [[ -f $DATA/$subject/rsfcPreprocessed/$subject.pm.cleanEPI.MNI.nii.gz ]] ; then
	cat <<EOF > run/${subject}_run_fwhm_acf.sh
#!/bin/bash
set -x
cd $DATA/${subject}/rsfcPreprocessed

# 3dcopy ${subject}.pm.mask.gm+orig. ${subject}.pm.mask.gm.nii

# applywarp \
# 	--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
# 	--in=${subject}.pm.mask.gm.nii.gz \
# 	--warp=../anat/${subject}.std.2.MNI.warpcoef.nii.gz \
# 	--out=${subject}.pm.mask.gm.MNI

#censor_file=$DATA/$subject/rsfcPreprocessed/${subject}.pm.censor.1D
#good_trs="\$( $scriptsDir/censor_to_good_list.r \$censor_file )"
#3dFWHMx -ACF -combine -detrend -mask ${subject}.mask.grey.MNI.nii.gz ${subject}.pm.cleanEPI.MNI.nii.gz"[\${good_trs}]" > ${subject}.cleanEPI.blur.est.1D
3dFWHMx -ACF -combine -detrend -mask ${subject}.mask.grey.MNI.nii.gz ${subject}.pm.cleanEPI.MNI.nii.gz > ${subject}.cleanEPI.blur.est.1D
 
EOF
	chmod u+x run/${subject}_run_fwhm_acf.sh
	echo "$scriptsDir/run/${subject}_run_fwhm_acf.sh" >> $taskFile
	## echo "cd $DATA/${subject}/rsfcPreprocessed; 3dFWHMx -ACF -combine -detrend -automask ${subject}.pm.cleanEPI.MNI.nii.gz > ${subject}.cleanEPI.blur.est.1D" >> $taskFile
    fi
done
exit
## jobname
#$ -N acf-est

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
exit
rm -f $ROOT/log/acf-est.log
sge_command="qsub -N acf-est -q all.q -j y -m n -V -wd $( pwd ) -o $ROOT/log/acf-est.log -t 1-$nTasks" 
echo $sge_command
( exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF
)
echo "Running qstat"
qstat
