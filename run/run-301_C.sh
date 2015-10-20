#!/bin/bash

## jobname
#$ -N aal-301_C

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

CPUS=1
JOBS="-jobs $CPUS"
OMP_NUM_THREADS=$CPUS

export JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT

cd /data/sanDiego/rsfcGraphAnalysis/scripts

#./01-preprocessAnatomy.sh    --subject=301_C
#./02-preprocessFunctional.sh --subject=301_C --drop=0 --fwhm=4.2 -c 0.3

./04-extractRsfcTimeseriesFromAalMasks.sh  --subject=301_C -l ../standard/aal_for_SPM8/fcseedlist3mm.txt 
## ./04-extractEstopTimeseriesFromAalMasks.sh --subject=301_C --seedlist ../standard/aal_for_SPM8/fcseedlist3mm.txt 

