#!/bin/bash

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/rsfcGraphAnalysis/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results
VBM_RESULTS=$DATA/vbm.subject.list.from.matthew.n111/stats/
MDD_STANDARD=$ROOT/standard/yeo17liberal/
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

task=rsfc

#ctrlSubjects="$( cat ../data/config/mdd.at.pilot.txt )"
#mddSubjects="$( cat ../data/config/mdd.nat.pilot.txt )"
subjects="$ctrlSubjects $mddSubjects"
#subjects="$( cat ../data/config/telomeresubjlist.txt )"
#subjects="$( cat testsubj.txt )"


#maskFile=$GROUP_RESULTS/estop_success-failures_leftmayg_mask+tlrc.BRIK.gz
#maskFile2=$GROUP_RESULTS/estop_success-failures_rightmayg_mask+tlrc.BRIK.gz

#maskFile=$GROUP_RESULTS/ttest-atANdNat.R_CMA+tlrc.BRIK.gz
#maskFile2=$GROUP_RESULTS/ttest_atAndNat.L_CMA+tlrc.BRIK.gz

#emotion=R_CMA.weight.3mm
#emotion2=L_CMA.weight.3mm

vbmFile=$VBM_RESULTS/GM_mod_merg_s3.nii.gz
maskFile=$MDD_STANDARD/yeo_17liberal_vbm.txt


for subject in $subjects ; do
    if [ ! -f $DATA/${subject}/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then

#roistatsFile=$GROUP_RESULTS/ohbm.fear-happy.${emotion}Failures.txt
#       cd $GROUP_RESULTS
#3dROIstats -nobriklab -mask $maskFile $DATA/$subject/functional/${subject}.estop.${emotion}Failures.%cs+tlrc.BRIK.gz > tempfile
#       sed '1D' tempfile >> $roistatsFile

roistatsFile2=$GROUP_RESULTS/k99pilot.${task}.${emotion}.txt
      cd $GROUP_RESULTS
3dROIstats -nobriklab -mask $maskFile $DATA/$subject/${task}/${emotion}/${emotion}.z-score+tlrc.BRIK.gz > tempfile
      sed '1D' tempfile >> $roistatsFile2

#roistatsFile3=$GROUP_RESULTS/ohbm.success-failures.${emotion}Success.txt
#cd $GROUP_RESULTS
#3dROIstats -nobriklab -mask $maskFile2 $DATA/$subject/functional/${subject}.estop.${emotion}Success.%cs+tlrc.BRIK.gz > tempfile
#sed '1D' tempfile >> $roistatsFile3

roistatsFile4=$GROUP_RESULTS/k99pilot.${task}.${emotion2}.txt
cd $GROUP_RESULTS
3dROIstats -nobriklab -mask $maskFile2 $DATA/$subject/${task}/${emotion2}/${emotion2}.z-score+tlrc.BRIK.gz > tempfile
sed '1D' tempfile >> $roistatsFile4

       else
           echo "**** $subject MOVED TOO MUCH ***"
      fi
done