#!/bin/bash

#set -x

trap exit SIGHUP SIGINT SIGTERM

taskDir=rsfcGraphAnalysis
task=ESTOP
## task=funcon

RSFC_ROOT=/data/sanDiego

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else
     # subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.nat.txt )"
     # subjects="106_C 112_C 118_C 144_C 149_C 158_C 161_C 300_C 304_C
     # 	      311_C 315_C 317_C 322_C 330_C 337_C 341_C 357_C 365_C
     # 	      367_C 389_C 397_C 403_C 410_C 415_C 111_C 117_C 120_C
     # 	      147_C 150_C 160_C 167_C 301_C 309_C 313_C 316_C 320_C
     # 	      323_C 336_C 339_C 348_C 364_C 366_C 380_C 392_C 401_C
     # 	      406_C 414_C 419_C"

    subjects="386_D 387_D 388_D 389_D 392_D 393_D 394_D 394_D2 397_D
    	      401_D 402_D 403_D 404_D 405_D 406_D 413_D 415_D 420_D 421_C 422_C
	      423_C 423_D 424_C 426_A"
    ## change task above too to match this list of subjects
##subjects=$( cd $RSFC_ROOT/$taskDir/data; \ls -1d *_[ABC] )
    ## subjects="$( cat /data/sanDiego/cESTOP/data/config/clean.estop.subjList.txt /data/sanDiego/cESTOP/data/config/clean.estop.mddfollowup.txt  )"
fi


date

if [[ ! -d ../log ]] ; then 
    mkdir ../log 
fi
	    
for subject in $subjects ; do
	
    ## if [[ -f $RSFC_ROOT/$taskDir/data/$subject/${subject}$task+orig.HEAD ]] ; then

	cat <<EOF > run/run-followup-$subject.sh
#!/bin/bash

## jobname
#$ -N followup-$subject

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

./01-preprocessAnatomy.sh    --subject=$subject

./02-makeNuisanceMasks.sh -s ${subject}

## ./02-preprocessFunctional.sh --subject=$subject --drop=0 --fwhm=4.2 -c 0.3

## ./04-extractRsfcTimeseriesFromAalMasks.sh  --subject=$subject -l ../standard/config/Harvard-Oxford_amygdala_seeds.txt

## ./04-extractEstopTimeseriesFromAalMasks.sh --subject=$subject --seedlist ../standard/aal_for_SPM8/fcseedlist3mm.txt 

EOF
	chmod +x  run/run-followup-$subject.sh
	echo $subject
	qsub -o ../log/followup-$subject.log \
	     run/run-followup-$subject.sh
    # else
    # 	echo "*** No such file: $RSFC_ROOT/$taskDir/data/$subject/${subject}$task+orig.HEAD"
    # fi
done
    
