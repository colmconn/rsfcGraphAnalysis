#!/bin/bash

set -x

if [[ 1 == 0 ]] ; then
    
## compute the union and intersection of the L_DLPFC from the between-group results
cd ../data/Group.results

3dcalc -a clorder.fwhm4.2.mddAndCtrl.L_BLA.weight.3mm+tlrc.HEAD -expr "equals(a, 2)" -prefix L_DLPFC.fwhm4.2.mddAndCtrl.L_BLA.weight.3mm
3dcalc -a clorder.fwhm4.2.mddAndCtrl.R_BLA.weight.3mm+tlrc.HEAD -expr "equals(a, 2)" -prefix L_DLPFC.fwhm4.2.mddAndCtrl.R_BLA.weight.3mm
3dcalc -a clorder.fwhm4.2.mddAndCtrl.L_SFA.weight.3mm+tlrc.HEAD -expr "equals(a, 6)" -prefix L_DLPFC.fwhm4.2.mddAndCtrl.L_SFA.weight.3mm

3dMean -mask_union -prefix mask.L_DLPFC.union \
       L_DLPFC.fwhm4.2.mddAndCtrl.L_BLA.weight.3mm+tlrc.HEAD \
       L_DLPFC.fwhm4.2.mddAndCtrl.R_BLA.weight.3mm+tlrc.HEAD \
       L_DLPFC.fwhm4.2.mddAndCtrl.L_SFA.weight.3mm+tlrc.HEAD

3dMean -mask_inter -prefix mask.L_DLPFC.intersection \
       L_DLPFC.fwhm4.2.mddAndCtrl.L_BLA.weight.3mm+tlrc.HEAD \
       L_DLPFC.fwhm4.2.mddAndCtrl.R_BLA.weight.3mm+tlrc.HEAD \
       L_DLPFC.fwhm4.2.mddAndCtrl.L_SFA.weight.3mm+tlrc.HEAD


## compute the union and intersection of the L_DLPFC from the CDRS-R regression results
cd ../Group.results.CDRS.t.score.scaled.diff.withAandC

3dcalc -a clorder.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       -expr "or(equals(a, 2), equals(a, 3))" \
       -prefix L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff

3dcalc -a clorder.regression.fwhm4.2.restingstate.mddOnly.R_SFA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       -expr "or(equals(a, 2), equals(a, 5))" \
       -prefix L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.R_SFA.weight.3mm.and.CDRS.t.score.scaled.diff

3dMean -mask_union -prefix mask.L_DLPFC.union \
       L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.R_SFA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD

3dMean -mask_inter -prefix mask.L_DLPFC.intersection \
       L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       L_DLPFC.regression.fwhm4.2.restingstate.mddOnly.R_SFA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD



## now compute the overlap of the between-group and CDRS-R regression L_DLPFC masks
3dcalc -a ../Group.results/mask.L_DLPFC.union+tlrc.HEAD \
       -b mask.L_DLPFC.union+tlrc.HEAD \
       -expr "a+2*b" \
       -prefix union.masks.L_DLPFC.from.between-group.and.CDRS-R.regression
3drefit -cmap INT_CMAP union.masks.L_DLPFC.from.between-group.and.CDRS-R.regression+tlrc.HEAD

3dcalc -a ../Group.results/mask.L_DLPFC.intersection+tlrc.HEAD \
       -b mask.L_DLPFC.intersection+tlrc.HEAD \
       -expr "a+2*b" \
       -prefix intersection.masks.L_DLPFC.from.between-group.and.CDRS-R.regression
3drefit -cmap INT_CMAP intersection.masks.L_DLPFC.from.between-group.and.CDRS-R.regression+tlrc.HEAD

3dclust -isomerge 1.01 0 union.masks.L_DLPFC.from.between-group.and.CDRS-R.regression+tlrc.HEAD        > clust.union.masks.L_DLPFC.from.between-group.and.CDRS-R.regression.txt
3dclust -isomerge 1.01 0 intersection.masks.L_DLPFC.from.between-group.and.CDRS-R.regression+tlrc.HEAD > clust.intersection.L_DLPFC.from.between-group.and.CDRS-R.regression.txt

\@DiceMetric ../Group.results/mask.L_DLPFC.union+tlrc.HEAD mask.L_DLPFC.union+tlrc.HEAD -save_match -save_diff > dice_metric.union.masks.L_DLPFC.from.between-group.and.CDRS-R.regression.txt

fi

cd ../data/Group.results.CDRS.t.score.scaled.diff.withAandC

## now compute the union of the R DLPFC ROIS, separate by positive or negative correlation
3dcalc -a clorder.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       -expr "or(equals(a, 1), equals(a, 9))" \
       -prefix R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff

3dcalc -a clorder.regression.fwhm4.2.restingstate.mddOnly.L_SFA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       -expr "equals(a, 2)" \
       -prefix R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_SFA.weight.3mm.and.CDRS.t.score.scaled.diff

3dcalc -a clorder.regression.fwhm4.2.restingstate.mddOnly.R_BLA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc.HEAD \
       -expr "equals(a, 5)" \
       -prefix R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.R_BLA.weight.3mm.and.CDRS.t.score.scaled.diff

3dcalc -a R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_CMA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc \
       -b R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.L_SFA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc \
       -c R_DLPFC.regression.fwhm4.2.restingstate.mddOnly.R_BLA.weight.3mm.and.CDRS.t.score.scaled.diff+tlrc \
       -expr "step(a + b) + 2*c" \
       -prefix R_DLPFC.CDRS.t.score.scaled.diff
3drefit -cmap INT_CMAP R_DLPFC.CDRS.t.score.scaled.diff+tlrc.HEAD

