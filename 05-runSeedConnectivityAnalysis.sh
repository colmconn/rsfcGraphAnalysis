#!/bin/bash

#set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else
    subjects=$( cd $DATA ;  ls -d *_[ABC] )
fi

taskFile=./runSeedConnectivity-TaskFile
cat /dev/null > $taskFile

for subject in $subjects ; do

    preprocessedRsfcDir=$DATA/$subject/rsfcPreprocessed

    if [[ -f ${preprocessedRsfcDir}/${subject}.pm.cleanEPI.MNI.nii.gz ]] ; then 

	# echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/Fox_seed_list.txt" >> ${taskFile}

	## echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/salience_seed_list.txt" >> ${taskFile}

	
	## echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/juelich_amygdala_seeds.txt" >> ${taskFile}
	
	## echo "./04-weightedSingleSubjectRsfc.sh -s $subject -l ../data/config/juelich_amygdala_seeds_weights.txt" >> ${taskFile}

	echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/juelich_whole_amygdala_seeds.txt" >> ${taskFile}
	
	## echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/Harvard-Oxford_amygdala_seeds.txt" >> ${taskFile}

	## the next oen was used fro the R01 preliminary data analysis
	## echo "./04-singleSubjectRsfc.sh -s $subject -l ../data/config/AAL_amygdala_seeds.txt" >> ${taskFile}


	## echo "./04-weightedSingleSubjectRsfcWithCleanedSeedTimeseries.sh -s $subject -l ../data/config/juelich_left_amygdala_seeds_weights.txt" >> ${taskFile}
	## echo "./04-weightedSingleSubjectRsfcWithCleanedSeedTimeseries.sh -s $subject -l ../data/config/juelich_right_amygdala_seeds_weights.txt" >> ${taskFile}		
    fi

done

## jobname
#$ -N seedConnectivity

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

rm -f $ROOT/log/runSeedConnectivity.log
sge_command="qsub -N seedConnectivity -q all.q -j y -m n -V -wd $( pwd ) -o $ROOT/log/runSeedConnectivity.log -t 1-$nTasks" 
echo $sge_command
exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF

echo "Running qstat"
qstat
