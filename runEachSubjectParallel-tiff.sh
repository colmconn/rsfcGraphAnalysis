#!/bin/bash

#set -x

trap exit SIGHUP SIGINT SIGTERM

taskDir=rsfcGraphAnalysis
task=funcon
#task=ESTOP

RSFC_ROOT=/data/sanDiego

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else
#subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/estop.realigned.subj.txt )"

subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/estop.subjList.txt )"
#subjects="$( cat /data/sanDiego/rsfcGraphAnalysis/data/config/rs114a.subj.txt )"
echo $subjects

    ## change task above too to match this list of subjects
## subjects=$( cd $RSFC_ROOT/$taskDir/data; \ls -1d *_[ABC] )
    ## subjects="$( cat /data/sanDiego/cESTOP/data/config/clean.estop.subjList.txt /data/sanDiego/cESTOP/data/config/clean.estop.mddfollowup.txt  )"
fi


date

if [[ ! -d ../log ]] ; then 
    mkdir ../log 
fi
	    
for subject in $subjects ; do
	
     if [[ -f $RSFC_ROOT/$taskDir/data/$subject/${subject}$task+orig.HEAD ]] ; then

	cat <<EOF > run/run-yeo7estop-$subject.sh
#!/bin/bash

## jobname
#$ -N yeo7-estop-$subject

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
JOBS="-jobs \$CPUS"
OMP_NUM_THREADS=\$CPUS

export JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT

cd $RSFC_ROOT/$taskDir/scripts

##./01-preprocessAnatomy.sh    --subject=$subject

#./02-preprocessFunctional.sh --subject=$subject --drop=0 --fwhm=4.2 -c 0.3

## ./04-extractRsfcTimeseriesFromAalMasks.sh  --subject=$subject -l ../standard/aal_for_SPM8/fcseedlist3mm.txt 
 ./04-extractRsfcTimeseriesFromAalMasks-tiff.sh --subject=$subject --seedlist /data/sanDiego/rsfcGraphAnalysis/standard/yeo7liberal/yeo_7liberal_split2.txt

EOF
	chmod +x  run/run-yeo7estop-$subject.sh
	echo $subject
	qsub -o ../log/$subject.yeo7estop.log \
	     run/run-yeo7estop-$subject.sh
     else
     	echo "*** No such file: $RSFC_ROOT/$taskDir/data/$subject/${subject}$task+orig.HEAD"
     fi
done
    
