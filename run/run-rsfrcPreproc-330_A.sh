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

## ./02-makeNuisanceMasks.sh -s 330_A --force
./02-makeNuisanceMasks.sh -s 330_A

cd /data/sanDiego/rsfcGraphAnalysis/data/330_A

[[ -d rsfcPreprocessed ]] && rm -fr rsfcPreprocessed

/home/colmconn/bin/afni_restproc.py \
	-zeropadfirst "-I 30" \
        -wmmask /data/sanDiego/rsfcGraphAnalysis/data/330_A/anat/white.mask+orig \
        -vmask /data/sanDiego/rsfcGraphAnalysis/data/330_A/anat/ventricles.mask+orig \
	-anat /data/sanDiego/rsfcGraphAnalysis/data/330_A/anat/330_A.anat_struc_brain.nii.gz \
	-epi  /data/sanDiego/rsfcGraphAnalysis/data/330_A/330_Afuncon+orig \
	-script 330_A.rsfc.preproc.tcsh \
	-dest /data/sanDiego/rsfcGraphAnalysis/data/330_A/rsfcPreprocessed \
	-prefix 330_A.pm \
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

cd rsfcPreprocessed
tcsh -f 330_A.rsfc.preproc.tcsh

if [[ -f 330_A.pm.cleanEPI+orig.HEAD ]] ; then 

	3dcopy 330_A.pm.cleanEPI+orig 330_A.pm.cleanEPI.nii

	echo "*** Warping cleanedEPI to MNI152 space"

	3dbucket -prefix 330_A.func.native.std.nii 330_A.pm.cleanEPI+orig\[0\]

	3dresample -master ../anat/330_A.anat_struc.std.nii.gz -prefix ./330_A.func.std.resampled.nii \
    		-inset 330_A.func.native.std.nii.gz

	flirt -in 330_A.func.native.std -ref 330_A.func.std.resampled \
	    -out 330_A.func2anat.flirt -omat 330_A.func2anat.flirt.mat

	applywarp \
		--ref=/data/sanDiego/rsfcGraphAnalysis/standard/MNI152_T1_3mm.nii.gz \
		--in=330_A.pm.cleanEPI.nii.gz \
		--warp=../anat/330_A.std.2.MNI.warpcoef \
		--premat=330_A.func2anat.flirt.mat \
		--out=330_A.pm.cleanEPI.MNI

	motionThresholdPrecentage=0.2
	threshold=20
	numberOfCensoredVolumes=$( cat 330_A.pm.censor.1D | gawk '{a+=(1-$0)}END{print a}' )
	totalNumberOfVolumes=$( cat 330_A.pm.censor.1D | wc -l )
	cutoff=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

	if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then 
	    echo "*** A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than $threshold % (n=$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_330_A_${threshold}percent.txt
	    echo "*** WARNING: 330_A will not be analysed due to having more than $threshold % of their volumes censored."
	fi

	3dFWHMx -combine -detrend -automask 330_A.pm.cleanEPI.MNI.nii.gz > 330_A.cleanEPI.blur.est.1D
	ln -sf ../anat/330_A.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD .
	ln -sf ../anat/330_A.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz .

fi 
