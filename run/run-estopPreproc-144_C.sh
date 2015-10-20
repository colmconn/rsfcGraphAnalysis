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

## ./02-makeNuisanceMasks.sh -s 144_C --force
./02-makeNuisanceMasks.sh -s 144_C

cd /data/sanDiego/rsfcGraphAnalysis/data/144_C

[[ -d estopPreprocessed ]] && rm -fr estopPreprocessed

/home/colmconn/bin/afni_restproc.py \
	-zeropadfirst "-I 30" \
        -wmmask /data/sanDiego/rsfcGraphAnalysis/data/144_C/anat/white.mask+orig \
        -vmask /data/sanDiego/rsfcGraphAnalysis/data/144_C/anat/ventricles.mask+orig \
	-anat /data/sanDiego/rsfcGraphAnalysis/data/144_C/anat/144_C.anat_struc_brain.nii.gz \
	-epi  /data/sanDiego/rsfcGraphAnalysis/data/144_C/144_CESTOP+orig \
	-script 144_C.estop.preproc.tcsh \
	-dest /data/sanDiego/rsfcGraphAnalysis/data/144_C/estopPreprocessed \
	-prefix 144_C.pm \
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

tcsh -f 144_C.estop.preproc.tcsh

if [[ -f 144_C.pm.cleanEPI+orig.HEAD ]] ; then 

	3dcopy 144_C.pm.cleanEPI+orig 144_C.pm.cleanEPI.nii

	echo "*** Warping cleanedEPI to MNI152 space"

	3dbucket -prefix 144_C.func.native.std.nii 144_C.pm.cleanEPI+orig\[0\]

	3dresample -master ../anat/144_C.anat_struc.std.nii.gz -prefix ./144_C.func.std.resampled.nii \
    		-inset 144_C.func.native.std.nii.gz

	flirt -in 144_C.func.native.std -ref 144_C.func.std.resampled \
	    -out 144_C.func2anat.flirt -omat 144_C.func2anat.flirt.mat

	applywarp \
		--ref=/data/sanDiego/rsfcGraphAnalysis/standard/MNI152_T1_3mm.nii.gz \
		--in=144_C.pm.cleanEPI.nii.gz \
		--warp=../anat/144_C.std.2.MNI.warpcoef \
		--premat=144_C.func2anat.flirt.mat \
		--out=144_C.pm.cleanEPI.MNI

	motionThresholdPrecentage=0.2
	threshold=20
	numberOfCensoredVolumes=$( cat 144_C.pm.censor.1D | gawk '{a+=(1-$0)}END{print a}' )
	totalNumberOfVolumes=$( cat 144_C.pm.censor.1D | wc -l )
	cutoff=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

	if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then 
	    echo "*** A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than $threshold % (n=$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_144_C_percent.txt
	    echo "*** WARNING: 144_C will not be analysed due to having more than $threshold % of their volumes censored."
	fi


	3dFWHMx -combine -detrend -automask 144_C.pm.cleanEPI.MNI.nii.gz > 144_C.cleanEPI.blur.est.1D
	ln -sf ../anat/144_C.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD .
	ln -sf ../anat/144_C.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz .

fi 
