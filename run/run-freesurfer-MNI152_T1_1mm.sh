#!/bin/bash

#$ -S /bin/bash

FREESURFER_HOME=/data/software/freesurfer/
SUBJECTS_DIR=/data/sanDiego/rsfcGraphAnalysis/data/
export FS_FREESURFERENV_NO_OUTPUT FREESURFER_HOME SUBJECTS_DIR
source /data/software/freesurfer//SetUpFreeSurfer.sh

recon-all -s MNI_152_T1_1mm -all

