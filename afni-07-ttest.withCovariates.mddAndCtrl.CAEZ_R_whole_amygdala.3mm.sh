#!/bin.bash

unset AFNI_COMPRESSOR 
export OMP_NUM_THREADS=40
mkdir -p /data/sanDiego/rsfcGraphAnalysis/data//Group.results.afni/new.ttest.CAEZ_R_whole_amygdala.3mm 
cd /data/sanDiego/rsfcGraphAnalysis/data//Group.results.afni/new.ttest.CAEZ_R_whole_amygdala.3mm 
3dttest++ -resid new.ttest.mddAndCtrl.CAEZ_R_whole_amygdala.3mm.covaried.resid \
-CLUSTSIM \
-prefix_clustsim resid.CStemp.CAEZ_R_whole_amygdala.3mm \
-mask /data/sanDiego/rsfcGraphAnalysis/standard//MNI_caez_N27_brain.3mm+tlrc.HEAD \
-prefix new.ttest.mddAndCtrl.CAEZ_R_whole_amygdala.3mm.covaried \
-center NONE \
-setA MDD 111_A /data/sanDiego/rsfcGraphAnalysis//data/111_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
112_A /data/sanDiego/rsfcGraphAnalysis//data/112_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
113_A /data/sanDiego/rsfcGraphAnalysis//data/113_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
118_A /data/sanDiego/rsfcGraphAnalysis//data/118_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
120_A /data/sanDiego/rsfcGraphAnalysis//data/120_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
133_A /data/sanDiego/rsfcGraphAnalysis//data/133_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
134_A /data/sanDiego/rsfcGraphAnalysis//data/134_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
136_A /data/sanDiego/rsfcGraphAnalysis//data/136_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
140_A /data/sanDiego/rsfcGraphAnalysis//data/140_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
144_A /data/sanDiego/rsfcGraphAnalysis//data/144_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
147_A /data/sanDiego/rsfcGraphAnalysis//data/147_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
154_A /data/sanDiego/rsfcGraphAnalysis//data/154_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
158_A /data/sanDiego/rsfcGraphAnalysis//data/158_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
160_A /data/sanDiego/rsfcGraphAnalysis//data/160_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
164_A /data/sanDiego/rsfcGraphAnalysis//data/164_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
167_A /data/sanDiego/rsfcGraphAnalysis//data/167_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
301_A /data/sanDiego/rsfcGraphAnalysis//data/301_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
304_A /data/sanDiego/rsfcGraphAnalysis//data/304_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
310_A /data/sanDiego/rsfcGraphAnalysis//data/310_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
330_A /data/sanDiego/rsfcGraphAnalysis//data/330_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
331_A /data/sanDiego/rsfcGraphAnalysis//data/331_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
332_A /data/sanDiego/rsfcGraphAnalysis//data/332_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
335_A /data/sanDiego/rsfcGraphAnalysis//data/335_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
336_A /data/sanDiego/rsfcGraphAnalysis//data/336_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
339_A /data/sanDiego/rsfcGraphAnalysis//data/339_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
342_A /data/sanDiego/rsfcGraphAnalysis//data/342_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
343_A /data/sanDiego/rsfcGraphAnalysis//data/343_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
351_A /data/sanDiego/rsfcGraphAnalysis//data/351_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
356_A /data/sanDiego/rsfcGraphAnalysis//data/356_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
358_A /data/sanDiego/rsfcGraphAnalysis//data/358_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
359_A /data/sanDiego/rsfcGraphAnalysis//data/359_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
360_A /data/sanDiego/rsfcGraphAnalysis//data/360_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
362_A /data/sanDiego/rsfcGraphAnalysis//data/362_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
364_A /data/sanDiego/rsfcGraphAnalysis//data/364_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
368_A /data/sanDiego/rsfcGraphAnalysis//data/368_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
371_A /data/sanDiego/rsfcGraphAnalysis//data/371_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
373_A /data/sanDiego/rsfcGraphAnalysis//data/373_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
378_A /data/sanDiego/rsfcGraphAnalysis//data/378_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
389_A /data/sanDiego/rsfcGraphAnalysis//data/389_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
393_A /data/sanDiego/rsfcGraphAnalysis//data/393_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
-setB NCL 108_A /data/sanDiego/rsfcGraphAnalysis//data/108_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
109_A /data/sanDiego/rsfcGraphAnalysis//data/109_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
116_A /data/sanDiego/rsfcGraphAnalysis//data/116_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
121_A /data/sanDiego/rsfcGraphAnalysis//data/121_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
122_A /data/sanDiego/rsfcGraphAnalysis//data/122_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
123_A /data/sanDiego/rsfcGraphAnalysis//data/123_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
126_A /data/sanDiego/rsfcGraphAnalysis//data/126_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
127_A /data/sanDiego/rsfcGraphAnalysis//data/127_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
131_A /data/sanDiego/rsfcGraphAnalysis//data/131_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
138_A /data/sanDiego/rsfcGraphAnalysis//data/138_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
139_A /data/sanDiego/rsfcGraphAnalysis//data/139_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
141_A /data/sanDiego/rsfcGraphAnalysis//data/141_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
142_A /data/sanDiego/rsfcGraphAnalysis//data/142_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
143_A /data/sanDiego/rsfcGraphAnalysis//data/143_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
145_A /data/sanDiego/rsfcGraphAnalysis//data/145_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
151_A /data/sanDiego/rsfcGraphAnalysis//data/151_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
152_A /data/sanDiego/rsfcGraphAnalysis//data/152_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
155_A /data/sanDiego/rsfcGraphAnalysis//data/155_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
157_A /data/sanDiego/rsfcGraphAnalysis//data/157_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
159_A /data/sanDiego/rsfcGraphAnalysis//data/159_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
162_A /data/sanDiego/rsfcGraphAnalysis//data/162_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
163_A /data/sanDiego/rsfcGraphAnalysis//data/163_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
165_A /data/sanDiego/rsfcGraphAnalysis//data/165_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
306_A /data/sanDiego/rsfcGraphAnalysis//data/306_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
307_A /data/sanDiego/rsfcGraphAnalysis//data/307_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
312_A /data/sanDiego/rsfcGraphAnalysis//data/312_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
321_A /data/sanDiego/rsfcGraphAnalysis//data/321_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
326_A /data/sanDiego/rsfcGraphAnalysis//data/326_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
328_A /data/sanDiego/rsfcGraphAnalysis//data/328_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
337_A /data/sanDiego/rsfcGraphAnalysis//data/337_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
354_A /data/sanDiego/rsfcGraphAnalysis//data/354_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
367_A /data/sanDiego/rsfcGraphAnalysis//data/367_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
374_A /data/sanDiego/rsfcGraphAnalysis//data/374_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
377_A /data/sanDiego/rsfcGraphAnalysis//data/377_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
380_A /data/sanDiego/rsfcGraphAnalysis//data/380_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
382_A /data/sanDiego/rsfcGraphAnalysis//data/382_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
395_A /data/sanDiego/rsfcGraphAnalysis//data/395_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
397_A /data/sanDiego/rsfcGraphAnalysis//data/397_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
400_A /data/sanDiego/rsfcGraphAnalysis//data/400_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
402_A /data/sanDiego/rsfcGraphAnalysis//data/402_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
406_A /data/sanDiego/rsfcGraphAnalysis//data/406_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
409_A /data/sanDiego/rsfcGraphAnalysis//data/409_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
414_A /data/sanDiego/rsfcGraphAnalysis//data/414_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
417_A /data/sanDiego/rsfcGraphAnalysis//data/417_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
423_A /data/sanDiego/rsfcGraphAnalysis//data/423_A/afniRsfc/CAEZ_R_whole_amygdala.3mm/CAEZ_R_whole_amygdala.3mm.z-score+tlrc.HEAD \
-covariates /data/sanDiego/rsfcGraphAnalysis/data//Group.data.afni/3dttest.covariates.mddAndCtrl.CAEZ_R_whole_amygdala.3mm.txt