#!/bin/bash

## jobname
#$ -N reg-RADS.Total.tscore

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
#JOBS="-jobs $CPUS"
#OMP_NUM_THREADS=$CPUS

export JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT

cd /data/sanDiego/rsfcGraphAnalysis/scripts

## ./rsfc.parallel.robust.regression.change.r  -c new.mdd.RADS.Total.tscore.diff.change.score.csv -v RADS.Total.tscore.diff -s juelich_amygdala_seeds_weights.txt

./rsfc.parallel.robust.regression.change.r  -c new.mdd.RADS.Total.tscore.scaled.diff.change.score.csv -v RADS.Total.tscore.scaled.diff -s juelich_amygdala_seeds_weights.txt

