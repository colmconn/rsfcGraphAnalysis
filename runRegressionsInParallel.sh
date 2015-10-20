#!/bin/bash

#set -x

trap exit SIGHUP SIGINT SIGTERM

taskDir=rsfcGraphAnalysis
task=funcon

RSFC_ROOT=/data/sanDiego

##regressionVariables="CGAS BDI.II.Total CDRS.t.score MASC.tscore CDI.Total RADS.Total.tscore"

regressionVariables="CDRS.t.score"

date

if [[ ! -d ../log ]] ; then 
    mkdir ../log 
fi
	    
for variable in $regressionVariables ; do
	
    cat <<EOF > run/run-${variable}-regression.sh
#!/bin/bash

## jobname
#$ -N reg-${variable}

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

AFNI_COMPRESSOR=GZIP
AFNI_DECONFLICT=OVERWRITE

#CPUS=1
#JOBS="-jobs \$CPUS"
#OMP_NUM_THREADS=\$CPUS

export JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT

cd $RSFC_ROOT/$taskDir/scripts

## ./rsfc.parallel.robust.regression.change.r  -e -c new.mdd.${variable}.diff.change.score.csv -v ${variable}.diff -s juelich_amygdala_seeds_weights.txt

./rsfc.parallel.robust.regression.change.r  -c new.mdd.${variable}.scaled.diff.change.score.csv -v ${variable}.scaled.diff -s juelich_amygdala_seeds_weights.txt

EOF
	chmod +x  run/run-${variable}-regression.sh
	echo $variable
	## qsub -pe smp 8 \
	rm -f ../log/${variable}-reversed-regression.log
	qsub -o ../log/${variable}-reversed-regression.log \
	     run/run-${variable}-regression.sh
done


    
#./rsfc.parallel.robust.regression.change.r  -c new.mdd.CDRS.t.score.rstandard.score.csv -v CDRS.t.score.rstandard -s juelich_amygdala_seeds_weights.txt
