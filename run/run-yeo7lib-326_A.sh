#!/bin/bash

## jobname
#$ -N yeo-326_A

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

##./01-preprocessAnatomy.sh    --subject=326_A

#./02-preprocessFunctional.sh --subject=326_A --drop=0 --fwhm=4.2 -c 0.3

## ./04-extractRsfcTimeseriesFromAalMasks.sh  --subject=326_A -l ../standard/aal_for_SPM8/fcseedlist3mm.txt 
 ./04-extractRsfcTimeseriesFromAalMasks-tiff.sh --subject=326_A --seedlist /data/sanDiego/rsfcGraphAnalysis/standard/yeo7liberal/yeo_7liberal_split.txt

