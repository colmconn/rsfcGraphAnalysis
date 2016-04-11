#!/bin/bash

#set -x

trap exit SIGHUP SIGINT SIGTERM

taskDir=rsfcGraphAnalysis
task=funcon

RSFC_ROOT=/data/sanDiego

##regressionVariables="CGAS BDI.II.Total CDRS.t.score MASC.tscore CDI.Total RADS.Total.tscore"

regressionVariables="CDRS.t.score"

## regressionVariables="CDRS.t.score.rstandard.short"



date

if [[ ! -d ../log ]] ; then 
    mkdir ../log 
fi
	    
for variable in $regressionVariables ; do
	
    cat <<EOF > run/run-${variable}-booted-regression.sh
#!/bin/bash

## jobname
#$ -N reg-${variable}

## queue
#$ -q all.q
_weights
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

AFNI_COMPRESSOR=GZIP
AFNI_DECONFLICT=OVERWRITE

#CPUS=1
#JOBS="-jobs \$CPUS"
#OMP_NUM_THREADS=\$CPUS

export JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT

cd $RSFC_ROOT/$taskDir/scripts

####################################################################################################
### Follow-up regressions

## ./rsfc.parallel.robust.regression.change.r  -e -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s juelich_amygdala_seeds_weights.txt

## ./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s juelich_amygdala_seeds_weights.txt

./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s left_sgacc_seeds.txt
./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s right_sgacc_seed.txt

####################################################################################################
### Baseline regression comands

##./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_left_whole_amygdala_seed.txt
##./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s juelich_right_whole_amygdala_seed.txt

./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s left_whole_sgacc_seed.txt
./rsfc.parallel.robust.regression.change.r  -b -r 100 -c new.mdd.${variable}.scores.csv   -v ${variable}       -s right_whole_sgacc_seed.txt

EOF
	chmod +x  run/run-${variable}-booted-regression.sh
	echo $variable
	## qsub -pe smp 8 \
	rm -f ../log/${variable}-booted-regression.log
#	qsub -o ../log/${variable}-booted-regression.log \
#	     run/run-${variable}-booted-regression.sh
done


    
#./rsfc.parallel.robust.regression.change.r  -c new.mdd.CDRS.t.score.rstandard.score.short.csv    -v CDRS.t.score.rstandard   -s juelich_whole_amygdala_seeds.txt
#./rsfc.parallel.robust.regression.change.r  -c new.mdd.CDRS.t.score.scaled.diff.change.score.csv -v CDRS.t.score.scaled.diff -s juelich_whole_amygdala_seeds.txt
#./rsfc.parallel.robust.regression.change.r  -c new.mdd.CDRS.t.score.score.withMedicated.csv      -v CDRS.t.score.scaled.diff -s juelich_whole_amygdala_seeds.txt
