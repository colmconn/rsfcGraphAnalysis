#!/bin/bash

FREESURFER_HOME=/data/software/freesurfer/
SUBJECTS_DIR=/data/sanDiego/rsfcGraphAnalysis/data/

export FREESURFER_HOME SUBJECTS_DIR

outputScriptName=run/run-freesurfer-MNI152_T1_1mm.sh


[[ ! -f $SUBJECTS_DIR/MNI_152_T1_1mm/mri/orig/001.mgz ]] && \
    mkdir $SUBJECTS_DIR/MNI_152_T1_1mm/mri/orig && \
    mri_convert -it nii -ot mgz -i $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz -o $SUBJECTS_DIR/MNI_152_T1_1mm/mri/orig/001.mgz

cat <<EOF > $outputScriptName
#!/bin/bash

#$ -S /bin/bash

FREESURFER_HOME=$FREESURFER_HOME
SUBJECTS_DIR=$SUBJECTS_DIR
export FS_FREESURFERENV_NO_OUTPUT FREESURFER_HOME SUBJECTS_DIR
source $FREESURFER_HOME/SetUpFreeSurfer.sh

recon-all -s MNI_152_T1_1mm -all

EOF


qsub -N MNI152-freesurfer -q all.q -j y -m n -V -wd $SUBJECTS_DIR -o ../log/MNI152-freesurfer.log $outputScriptName
qstat

