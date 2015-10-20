#!/bin/bash

set -x 

#$ -S /bin/bash

export PYTHONPATH=/data/software/afni

## use the newer faster despiking method. comment this out to get the
## old one back
export AFNI_3dDespike_NEW=YES

## only use a single thread since we're going to run so many subjects
## in parallel
export OMP_NUM_THREADS=1

## ./02-makeNuisanceMasks.sh -s 389_A --force
./02-makeNuisanceMasks-tiff.sh -s 389_A

cd /data/sanDiego/rsfcGraphAnalysis/data/389_A

[[ -d estopPreprocessed ]] && rm -fr estopPreprocessed

/data/sanDiego/rsfcGraphAnalysis/scripts/afni_restproc.py \
	-zeropadfirst "-I 30" \
        -wmmask /data/sanDiego/rsfcGraphAnalysis/data/389_A/anat/white.mask+orig \
        -vmask /data/sanDiego/rsfcGraphAnalysis/data/389_A/anat/ventricles.mask+orig \
	-anat /data/sanDiego/rsfcGraphAnalysis/data/389_A/anat/389_A.anat_struc_brain.nii.gz \
	-epi  /data/sanDiego/rsfcGraphAnalysis/data/389_A/389_AESTOP+orig \
	-script 389_A.estop.preproc.tcsh \
	-dest /data/sanDiego/rsfcGraphAnalysis/data/389_A/estopPreprocessed \
	-prefix 389_A.pm \
	-despike on \
	-anat_has_skull no \
	-tsnr \
	-epi2anat \
	-align on \
	-trcut 4 \
	-dvarscensor \
	-episize 3 \
	-dreg \
	-smooth on \
	-smoothrad 4.2 \
	-smoothfirst \
	-bandpass \
	-bpassregs \
	-polort 2 \
	-censorleft 1 \
	-censorright 2 \
	-fdlimit 0.5 \
	-dvarslimit 5 \
	-globalnorm \
	-corrmap \
	-cleanOnly \
	-exec off

cd estopPreprocessed
tcsh -f 389_A.estop.preproc.tcsh

if [[ -f 389_A.pm.cleanEPI+orig.HEAD ]] ; then 

	3dcopy 389_A.pm.cleanEPI+orig 389_A.pm.cleanEPI.nii

	echo "*** Warping cleanedEPI to MNI152 space"

	3dbucket -prefix 389_A.func.native.std.nii 389_A.pm.cleanEPI+orig\[0\]

	3dresample -master ../anat/389_A.anat_struc.std.nii.gz -prefix ./389_A.func.std.resampled.nii \
    		-inset 389_A.func.native.std.nii.gz

	flirt -in 389_A.func.native.std -ref 389_A.func.std.resampled \
	    -out 389_A.func2anat.flirt -omat 389_A.func2anat.flirt.mat

	applywarp \
		--ref=/data/sanDiego/rsfcGraphAnalysis/standard/MNI152_T1_3mm.nii.gz \
		--in=389_A.pm.cleanEPI.nii.gz \
		--warp=../anat/389_A.std.2.MNI.warpcoef \
		--premat=389_A.func2anat.flirt.mat \
		--out=389_A.pm.cleanEPI.MNI

	motionThresholdPrecentage=0.2
	threshold=20
	numberOfCensoredVolumes=$( cat 389_A.pm.censor.1D | gawk '{a+=(1-$0)}END{print a}' )
	totalNumberOfVolumes=$( cat 389_A.pm.censor.1D | wc -l )
	cutoff=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

	if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then 
	    echo "*** A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than $threshold % (n=$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_389_A_${threshold}percent.txt
	    echo "*** WARNING: 389_A will not be analysed due to having more than $threshold % of their volumes censored."
	fi


	3dFWHMx -combine -detrend -automask 389_A.pm.cleanEPI.MNI.nii.gz > 389_A.cleanEPI.blur.est.1D
	ln -sf ../anat/389_A.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD .
	ln -sf ../anat/389_A.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz .

fi 
