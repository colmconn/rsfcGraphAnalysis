#!/bin/bash

## jobname
#$ -N yeo-329_A

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

##./01-preprocessAnatomy.sh    --subject=329_A

#./02-preprocessFunctional.sh --subject=329_A --drop=0 --fwhm=4.2 -c 0.3

## ./04-extractRsfcTimeseriesFromAalMasks.sh  --subject=329_A -l ../standard/aal_for_SPM8/fcseedlist3mm.txt 
 ./04-extractRsfcTimeseriesFromAalMasks-tiff.sh --subject=329_A --seedlist ../standard/yeo_rsns_17tight_separated/yeo_rsn_17tight.txt

