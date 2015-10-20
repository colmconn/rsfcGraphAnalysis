#!/bin/bash

set -x 

##FREESURFER_HOME=/data/software/freesurfer/
##SUBJECTS_DIR=/data/sanDiego/freesurferAnalysis/data

DATA=$( readlink -e ../data )

#subjects=$( cd $SUBJECTS_DIR ;  ls -d *_A )

subjects="105_A"

[[ -d run ]] || mkdir run

for subject in $subjects ; do
    echo "####################################################################################################"
    echo "Generating script for subject $subject"

    if  [[ -f /data/sanDiego/restingstate/data/$subject/${subject}funcon+orig.HEAD ]] && \
	[[ -f /data/sanDiego/restingstate/data/$subject/${subject}funcon+orig.BRIK.gz ]]  ; then

	rsfcFile=/data/sanDiego/restingstate/data/$subject/${subject}funcon+orig

    elif [[ -f /data/sanDiego/${subject}BRIKS/${subject}funcon+orig.HEAD ]] && \
	[[ -f /data/sanDiego/${subject}BRIKS/${subject}funcon+orig.BRIK.gz ]]  ; then

	rsfcFile=/data/sanDiego/${subject}BRIKS/${subject}funcon+orig

    else 
	echo "*** Can not find resting state EPI file for ${subject}. Skipping."
	break
    fi

if  [[ -f /data/sanDiego/restingstate/data/$subject/${subject}+orig.HEAD ]] && \
	[[ -f /data/sanDiego/restingstate/data/$subject/${subject}+orig.BRIK.gz ]]  ; then

	anatFile=/data/sanDiego/restingstate/data/$subject/${subject}+orig

    elif [[ -f /data/sanDiego/${subject}BRIKS/${subject}funcon+orig.HEAD ]] && \
	[[ -f /data/sanDiego/${subject}BRIKS/${subject}funcon+orig.BRIK.gz ]]  ; then

	anatFile=/data/sanDiego/${subject}BRIKS/${subject}+orig

    else 
	echo "*** Can not find anat T1file for ${subject}. Skipping."
	break
    fi

cat <<EOF > run/run-${subject}.sh
#!/bin/bash

set -x 

#$ -S /bin/bash

##SUBJECTS_DIR=$SUBJECTS_DIR
DATA=$DATA

subject=$subject
anatFile=$anatFile
rsfcFile=$rsfcFile

[[ -d $DATA/$subject ]] && rm -rf $DATA/$subject


mkdir $DATA/$subject
cd $DATA/$subject

ln -sf \${anatFile}.HEAD ./
ln -sf \${anatFile}.BRIK.gz ./

ln -sf \${rsfcFile}.HEAD ./
ln -sf \${rsfcFile}.BRIK.gz ./

## convert absolute to relative symlinks
symlinks -c ./

anatFile=${anatFile##*/}
rsfcFile=${rsfcFile##*/}

cd $DATA/$subject

3dUnifize -prefix \${anatFile%%+*}.anat.unifized \${anatFile}
3dSkullStrip -input \${anatFile%%+*}.anat.unifized+orig -prefix \${anatFile%%+*}.anat.noskull 
3dSeg    -anat \${anatFile%%+*}.anat.unifized+orig  -mask AUTO

# afni_restproc.py \\
# 	-anat \${anatFile%%+*}.anat.noskull+orig \\
# 	-epi  \$rsfcFile \\
# 	-script \${subject}.rsfc.preproc.tcsh \\
# 	-dest \$DATA/\$subject \\
# 	-prefix $subject.pm \\
# 	-despike on \\
# 	-anat_has_skull no \\
# 	-tsnr \\
# 	-epi2anat \\
# 	-align on \\
# 	-trcut 4 \\
# 	-dvarscensor \\
# 	-episize 3 \\
# 	-dreg \\
# 	-smooth on \\
# 	-smoothrad 4.2 \\
# 	-smoothfirst \\
# 	-bandpass \\
# 	-bpassregs \\
# 	-polort 2 \\
# 	-censorleft 1 \\
# 	-censorright 2 \\
# 	-fdlimit 0.5 \\
# 	-dvarslimit 5 \\
# 	-globalnorm \\
# 	-corrmap \\
# 	-exec off

echo "( cd $DATA/$subject ; tcsh -f ${subject}.rsfc.preproc.tcsh )"
echo "( cd $DATA/$subject/tmp; gzip -9 *.nii )"
 
EOF

chmod +x run/run-${subject}.sh
#echo "qsub -N rsfc-$subject -q all.q -j y -m n -V -wd $DATA -o rsfc-${subject}.log"

done



## version of preprocessing command for use with freesurfer processed data 

## smoothing of the fMRI data is turned off as this will be done in
## the freesurfer pipeline see 1. Yeo BTT, Krienen FM, Sepulcre J,
## Sabuncu MR, Lashkari D, et al. (2011) The organization of the human
## cerebral cortex estimated by intrinsic functional connectivity. J
## Neurophysiol 106: 1125–1165. doi:10.1152/jn.00338.2011.

# afni_restproc.py \\
# 	-aseg $SUBJECTS_DIR/$subject/mri/aseg.mgz \\
# 	-anat $SUBJECTS_DIR/$subject/mri/brainmask.mgz \\
# 	-epi $rsfcFile \\
# 	-script ${subject}.rsfc.preproc.tcsh \\
# 	-dest $DATA/$subject \\
# 	-prefix $subject.pm \\
# 	-despike on \\
# 	-anat_has_skull no \\
# 	-tsnr \\
# 	-epi2anat \\
# 	-align on \\
# 	-trcut 4 \\
# 	-dvarscensor \\
# 	-episize 3 \\
# 	-dreg \\
# 	-smooth off \\
# 	-bandpass \\
# 	-bpassregs \\
# 	-polort 2 \\
# 	-censorleft 1 \\
# 	-censorright 2 \\
# 	-fdlimit 0.5 \\
# 	-dvarslimit 5 \\
# 	-globalnorm \\
# 	-corrmap \\
# 	-exec off
