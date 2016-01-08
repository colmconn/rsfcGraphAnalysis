#!/bin/bash

#set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard

taskFile=./multipleParallBootedRegressions-TaskFile
cat /dev/null > $taskFile

#regressionVariables="CDRS.t.score.scaled.diff" ## Predictive regressions
## regressionVariables="CDRS.t.score" ## Baseline regressions
# regressionVariables="CDRS.t.score.C"

#regressionVariables="CDRS.t.score.rstandard"

# regressionVariables="BDI.II.Total"
regressionVariables="MASC.tscore"
for variable in $regressionVariables ; do

    ## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt"  >> ${taskFile}
    ## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt" >> ${taskFile}

    echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.timepoint.a.score.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt"  >> ${taskFile}
    echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.timepoint.a.score.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt" >> ${taskFile}

done

## jobname
#$ -N multipleParallBootedRegressions

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

rm -f $ROOT/log/multipleParallBootedRegressions.log
sge_command="qsub -N multipleParallBootedRegressions -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
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
