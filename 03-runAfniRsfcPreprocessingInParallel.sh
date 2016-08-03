#!/bin/bash

#set -x 

trap exit SIGHUP SIGINT SIGTERM

studyName=rsfcGraphAnalysis

ROOT=/data/sanDiego/$studyName
DATA=$ROOT/data
RAW_DATA=$DATA/raw
PROCESSED_DATA=$DATA
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts
SUBJECTS_DIR=$PROCESSED_DATA

# if [[ "$#" -gt 0 ]] ; then
#     subjects="$@"
# else
#     subjects=$( cd $DATA ;  ls -d *_[ABC] )
# fi

subjects="107_A 109_A 137_A 151_A 153_A 154_A 164_A 300_A 304_A 325_A 332_A 342_A 343_A 345_A 348_A 351_A 356_A 360_A 366_A 378_A 376_A 380_A 395_A"

[[ -d run ]] || mkdir run

for subject in $subjects ; do
    echo "####################################################################################################"
    echo "*** Generating script for subject $subject"

    if  [[ ! -f $DATA/$subject/${subject}funcon+orig.HEAD ]] && \
	[[ ! -f $DATA/$subject/${subject}funcon+orig.BRIK.gz ]]  ; then

	echo "*** Can not find resting state EPI file for ${subject}. Skipping."
	continue
    else
	epiFile=$DATA/$subject/${subject}funcon+orig.HEAD
    fi

    if  [[ ! -f ${SUBJECTS_DIR}/$subject/$subject.anat.nii.gz ]] ; then #  && \
	# [[ ! -f ${SUBJECTS_DIR}/$subject/$subject.anat+orig.BRIK.gz ]]  ; then

	echo "*** Can not find anatomy file for subject ${subject}. Skipping."
	continue
    else
	anatFile=${SUBJECTS_DIR}/$subject/$subject.anat.nii.gz
    fi

    outputScriptName=run/run-afniRsfcPreproc-${subject}.sh

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

cd $SUBJECTS_DIR/$subject

preprocessingScript=${subject}.afniRsfcPreprocess.csh
rm -f \${preprocessingScript}

outputDir=afniRsfcPreprocessed
rm -fr \${outputDir}

motionThreshold=0.2
outlierThreshold=0.1

## 	     -tlrc_NL_warp							\\
## the regress_ROI WMe is not needed when using anaticor to perform WM signal removal
##	     -regress_ROI WMe							\\
##           -align_opts_aea "-giant_move"                                      \\

afni_proc.py -subj_id ${subject}						\\
             -script \${preprocessingScript}					\\
	     -out_dir \${outputDir}						\\
	     -blocks despike tshift align tlrc volreg mask blur regress		\\
	     -copy_anat $anatFile                                               \\
	     -dsets $epiFile                                                    \\
	     -tcat_remove_first_trs 3						\\
	     -tlrc_base MNI_caez_N27+tlrc					\\
	     -volreg_align_to MIN_OUTLIER					\\
	     -volreg_tlrc_warp							\\
	     -blur_size 8                                                       \\
	     -blur_to_fwhm  							\\
	     -blur_opts_B2FW "-ACF -rate 0.2 -temper"                           \\
	     -mask_apply group							\\
	     -mask_segment_anat yes						\\
	     -mask_segment_erode yes						\\
	     -regress_anaticor							\\
	     -regress_bandpass 0.01 0.1						\\
	     -regress_apply_mot_types demean deriv				\\
             -regress_censor_motion \$motionThreshold              		\\
	     -regress_censor_outliers \$outlierThreshold                 	\\
	     -regress_run_clustsim no						\\
	     -regress_est_blur_errts

if [[ -f \${preprocessingScript} ]] ; then 
   tcsh -xef \${preprocessingScript}

    cd \${outputDir}
    xmat_regress=X.xmat.1D 

    if [[ -f \$xmat_regress ]] ; then 

	motionThresholdFraction=0.2
	threshold=20
        fractionOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts frac_cen )
        numberOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_cen )
        totalNumberOfVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_no_cen )

	## cutoff=\$( echo "scale=0; \$motionThresholdFraction*\$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )
        ## rounding method from http://www.alecjacobson.com/weblog/?p=256
        cutoff=\$( echo "(\$(echo "scale=0;\$motionThresholdFraction*\$totalNumberOfVolumes" | bc)+0.5)/1" | bc )
	if [[ \$numberOfCensoredVolumes -gt \$cutoff ]] ; then 
	    echo "*** A total of \$numberOfCensoredVolumes of \$totalNumberOfVolumes we censored which is greater than \$motionThresholdFraction  (n=\$cutoff) of all total volumes of this subject" > 00_DO_NOT_ANALYSE_${subject}_\${threshold}percent.txt
	    echo "*** WARNING: $subject will not be analysed due to having more than \$threshold % of their volumes censored."
	fi

    	trs=\$( 1d_tool.py -infile \$xmat_regress -show_trs_uncensored encoded   \\
                          -show_trs_run 01 )
    	if [[ \$trs != "" ]] ; then  
	    3dFWHMx -ACF -detrend -mask mask_group+tlrc                      \\
        	errts.$subject.anaticor+tlrc"[\$trs]" > acf.blur.errts.1D
	fi 
    fi
else
    echo "*** No such file \${preprocessingScript}"
    echo "*** Cannot continue"
    exit 1
fi	

EOF

    chmod +x $outputScriptName
    rm -f ../log/$subject-rsfc-afniPreproc.log
    qsub -N rsfc-$subject -q all.q -j y -m n -V -wd $( pwd )  -o ../log/$subject-rsfc-afniPreproc.log $outputScriptName

done

qstat

# freeview -v \
# mri/T1.mgz \
# mri/wm.mgz \
# mri/brainmask.mgz \
# mri/aparc.a2009s+aseg.mgz:colormap=lut:opacity=0.2 \
# -f surf/lh.white:edgecolor=blue \
# surf/lh.pial:edgecolor=red \
# surf/rh.white:edgecolor=blue \
# surf/rh.pial:edgecolor=red
