#!/bin/bash

##set -x

ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
task=ESTOP
scriptDir=$ROOT/scripts

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else
    # subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.subjectList.txt )"
    
    #   subjects="$( cat ../data/config/mddfollowup.txt )"
    #subjects=$( cd $DATA ;  ls -d *_[AC] )
    
    ##subjects="118_C 120_C 149_C 150_C 301_C 311_A 315_C 364_A 364_C"
    ## subjects="324_B"
    
    ##     subjects="$( cat /data/sanDiego/cESTOP/data/config/clean.estop.mddfollowup.txt /data/sanDiego/cESTOP/data/config/clean.estop.subjList.txt  )"
    subjects="$( cat /data/sanDiego/cESTOP/data/config/estop.subjList.txt  )"

fi



[[ -d run ]] || mkdir run

for subject in $subjects ; do
    echo "####################################################################################################"
    echo "Generating script for subject $subject"

    if  [[ -f $DATA/$subject/${subject}${task}+orig.HEAD ]] && \
	[[ -f $DATA/$subject/${subject}${task}+orig.BRIK.gz ]]  ; then

	estopPrefix=$DATA/$subject/${subject}${task}+orig
    else 
	echo "*** Can not find EPI file for ${subject}. Skipping."
	continue
    fi

    if  [[ -f $DATA/$subject/anat/${subject}.anat_struc_brain.nii.gz ]] ; then
    
	skullStrippedAnatPrefix=$DATA/$subject/anat/${subject}.anat_struc_brain.nii.gz
	wmmask=$DATA/$subject/anat/white.mask+orig
	vmask=$DATA/$subject/anat/ventricles.mask+orig
    else 
	echo "*** Can not find anat T1 anatomy file for ${subject}. Skipping."
	continue
    fi

    outputScriptName=run/run-estopPreproc-${subject}.sh

    cat <<EOF > $outputScriptName
#!/bin/bash

set -x 

#$ -S /bin/bash

export PYTHONPATH=$AFNI_R_DIR

## use the newer faster despiking method. comment this out to get the
## old one back
export AFNI_3dDespike_NEW=YES

## only use a single thread since we're going to run so many subjects
## in parallel
export OMP_NUM_THREADS=1

## ./02-makeNuisanceMasks.sh -s ${subject} --force
./02-makeNuisanceMasks-tiff.sh -s ${subject}

cd $DATA/$subject

[[ -d estopPreprocessed ]] && rm -fr estopPreprocessed

$scriptDir/afni_restproc.py \\
	-zeropadfirst "-I 30" \\
        -wmmask $wmmask \\
        -vmask $vmask \\
	-anat ${skullStrippedAnatPrefix} \\
	-epi  $estopPrefix \\
	-script ${subject}.estop.preproc.tcsh \\
	-dest $DATA/$subject/estopPreprocessed \\
	-prefix $subject.pm \\
	-despike on \\
	-anat_has_skull no \\
	-tsnr \\
	-epi2anat \\
	-align on \\
	-trcut 4 \\
	-dvarscensor \\
	-episize 3 \\
	-dreg \\
	-smooth on \\
	-smoothrad 4.2 \\
	-smoothfirst \\
	-bandpass \\
	-bpassregs \\
	-polort 2 \\
	-censorleft 1 \\
	-censorright 2 \\
	-fdlimit 0.5 \\
	-dvarslimit 5 \\
	-globalnorm \\
	-corrmap \\
	-cleanOnly \\
	-exec off

cd estopPreprocessed
EOF
    
    if [[ ${subject} == 364_A ]] ; then
	cat <<EOF >> $outputScriptName

cat ${subject}.estop.preproc.tcsh | sed -e "s/align_epi_anat.py/align_epi_anat.py -giant_move/" >  ${subject}.estop.preproc.tcsh.new
mv -f ${subject}.estop.preproc.tcsh.new ${subject}.estop.preproc.tcsh

EOF
    fi

    cat <<EOF >> $outputScriptName
tcsh -f ${subject}.estop.preproc.tcsh

if [[ -f ${subject}.pm.cleanEPI+orig.HEAD ]] ; then 

	3dcopy ${subject}.pm.cleanEPI+orig ${subject}.pm.cleanEPI.nii

	echo "*** Warping cleanedEPI to MNI152 space"

	3dbucket -prefix ${subject}.func.native.std.nii ${subject}.pm.cleanEPI+orig\\[0\\]

	3dresample -master ../anat/${subject}.anat_struc.std.nii.gz -prefix ./$subject.func.std.resampled.nii \\
    		-inset ${subject}.func.native.std.nii.gz

	flirt -in ${subject}.func.native.std -ref ${subject}.func.std.resampled \\
	    -out ${subject}.func2anat.flirt -omat ${subject}.func2anat.flirt.mat

	applywarp \\
		--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \\
		--in=${subject}.pm.cleanEPI.nii.gz \\
		--warp=../anat/${subject}.std.2.MNI.warpcoef \\
		--premat=${subject}.func2anat.flirt.mat \\
		--out=${subject}.pm.cleanEPI.MNI

	motionThresholdPrecentage=0.2
	threshold=20
	numberOfCensoredVolumes=\$( cat ${subject}.pm.censor.1D | gawk '{a+=(1-\$0)}END{print a}' )
	totalNumberOfVolumes=\$( cat ${subject}.pm.censor.1D | wc -l )
	cutoff=\$( echo "scale=0; \$motionThresholdPrecentage*\$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

	if [[ \$numberOfCensoredVolumes -gt \$cutoff ]] ; then 
	    echo "*** A total of \$numberOfCensoredVolumes of \$totalNumberOfVolumes we censored which is greater than \$threshold % (n=\$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_${subject}_\${threshold}percent.txt
	    echo "*** WARNING: $subject will not be analysed due to having more than \$threshold % of their volumes censored."
	fi


	3dFWHMx -combine -detrend -automask ${subject}.pm.cleanEPI.MNI.nii.gz > ${subject}.cleanEPI.blur.est.1D
	ln -sf ../anat/${subject}.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD .
	ln -sf ../anat/${subject}.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz .

fi 
EOF

    chmod +x $outputScriptName
    echo -n "Submitting job..."
    rm -f ../log/$subject-estop-preproc.log
    qsub -N estop-$subject -q all.q -j y -m n -V -wd $( pwd )  -o ../log/$subject-estop-preproc.log $outputScriptName

done

