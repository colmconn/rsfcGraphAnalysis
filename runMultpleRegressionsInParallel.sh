#!/bin/bash

#set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard

taskFile=./multipleParallelBootedRegressions-TaskFile
cat /dev/null > $taskFile

#regressionVariables="CDRS.t.score.scaled.diff" ## Predictive regressions
regressionVariables="CDRS.t.score" ## Baseline regressions
# regressionVariables="CDRS.t.score.C"

#regressionVariables="CDRS.t.score.rstandard"

# regressionVariables="BDI.II.Total"
## regressionVariables="MASC.tscore"
for variable in $regressionVariables ; do

    ## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt"  >> ${taskFile}
    ## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt" >> ${taskFile}

##    echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.timepoint.a.score.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt"  >> ${taskFile}
##    echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.timepoint.a.score.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt" >> ${taskFile}


####################################################################################################
### Follow-up regressions

## ./rsfc.parallel.robust.regression.change.r  -e -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s juelich_amygdala_seeds_weights.txt >> ${taskFile}

## ./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s juelich_amygdala_seeds_weights.txt >> ${taskFile}


    ## Jan 17, 2017: The following two lines were used to generate the
    ## regression results reported in Connolly, C. G. et
    ## al. Resting-state functional connectivity of the amygdala and
    ## longitudinal changes in depression severity in adolescent
    ## depression. J Affect Disord 207, 86–94 (2017).
    
## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s left_whole_sgacc_seed.txt" >> ${taskFile}
## echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s right_whole_sgacc_seed.txt" >> ${taskFile}

## extra seed runs go here

# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s     juelich_whole_amygdala_seeds.txt" >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s     short_ACC_seed_list.txt" >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s     hippocampus_ventral_striatum_seeds.txt" >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s     followup-dlpfc-ins-IP-MPFC-seeds.txt" >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s     Fox-Goldapple-seeds.txt" >> ${taskFile}

## these are used for the predictive analysis that will got in Tony's R01 submission
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     juelich_whole_amygdala_seeds.txt"       >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     short_ACC_seed_list.txt"                >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     hippocampus_ventral_striatum_seeds.txt" >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     followup-dlpfc-ins-IP-MPFC-seeds.txt"   >> ${taskFile}
# echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     Fox-Goldapple-seeds.txt"                >> ${taskFile}

echo "./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s     miller-dmn.txt"                >> ${taskFile}


####################################################################################################
### Baseline regression comands

##./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt >> ${taskFile}
##./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt >> ${taskFile}

#echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s left_whole_sgacc_seed.txt" >> ${taskFile}
#echo "./rsfc.parallel.robust.regression.change.r  -b -r 1000 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s right_whole_sgacc_seed.txt" >> ${taskFile}


    
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

rm -f $ROOT/log/multipleParallelBootedRegressions.log
sge_command="qsub -N multipleParallelBootedRegressions -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
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
