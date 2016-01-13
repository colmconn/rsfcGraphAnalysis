#!/bin/bash

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

# set -x

cd ../data/seeds

[[ ! -d sgACC ]] && mkdir sgACC

cd sgACC

3dcalc -a $FSLDIR/data/atlases/Talairach/Talairach-labels-1mm.nii.gz -expr "equals(a, 381)" -prefix R.sgACC

# derived from the number of voxels in the ROI as extracted from the
# 1mm Tailarach atlas included with FSL
nvoxels=388
3dcalc -a R.sgACC+tlrc -expr "a+0.01" -prefix R.sgACC.plusSmallValue -datum float
3dmerge -1blur_fwhm 4.2 -prefix R.sgACC.plusSmallValue.blurred R.sgACC.plusSmallValue+tlrc
3dclust -1Dformat -savemask R.sgACC.plusSmallValue.clustered -nosum -1dindex 0 -1tindex 0 -2thresh -0.029 0.029 -dxyz=1 1.01 $nvoxels R.sgACC.plusSmallValue.blurred+tlrc > clust.R.sgACC.txt
@clip_volume -input R.sgACC.plusSmallValue.clustered+tlrc -left 0 

mv -f R.sgACC.plusSmallValue.clustered_clp+tlrc.HEAD R.sgACC.mask+tlrc.HEAD
mv -f R.sgACC.plusSmallValue.clustered_clp+tlrc.BRIK.gz R.sgACC.mask+tlrc.BRIK.gz

3dresample -rmode NN -master ../../../standard/MNI152_T1_3mm_brain.nii.gz -inset R.sgACC.mask+tlrc.HEAD -prefix R.sgACC.mask.3mm



####################################################################################################
3dcalc -a $FSLDIR/data/atlases/Talairach/Talairach-labels-1mm.nii.gz -expr "equals(a, 380)" -prefix L.sgACC

# derived from the number of voxels in the ROI as extracted from the
# 1mm Tailarach atlas included with FSL
nvoxels=441
3dcalc -a L.sgACC+tlrc -expr "a+0.01" -prefix L.sgACC.plusSmallValue -datum float
3dmerge -1blur_fwhm 4.2 -prefix L.sgACC.plusSmallValue.blurred L.sgACC.plusSmallValue+tlrc
3dclust -1Dformat -savemask L.sgACC.plusSmallValue.clustered -nosum -1dindex 0 -1tindex 0 -2thresh -0.029 0.029 -dxyz=1 1.01 $nvoxels L.sgACC.plusSmallValue.blurred+tlrc > clust.L.sgACC.txt
@clip_volume -input L.sgACC.plusSmallValue.clustered+tlrc -right 0 

mv -f L.sgACC.plusSmallValue.clustered_clp+tlrc.HEAD L.sgACC.mask+tlrc.HEAD
mv -f L.sgACC.plusSmallValue.clustered_clp+tlrc.BRIK.gz L.sgACC.mask+tlrc.BRIK.gz


3dresample -rmode NN -master ../../../standard/MNI152_T1_3mm_brain.nii.gz -inset L.sgACC.mask+tlrc.HEAD -prefix L.sgACC.mask.3mm


## create the combined mask

3dcalc -a R.sgACC.mask+tlrc     -b L.sgACC.mask+tlrc     -expr "a+2*b" -prefix sgACC.mask
3dcalc -a R.sgACC.mask.3mm+tlrc -b L.sgACC.mask.3mm+tlrc -expr "a+2*b" -prefix sgACC.mask.3mm

3dclust -isomerge 0 0  sgACC.mask.3mm+tlrc > clust.sgACC.3mm.txt

ln -sf ../MNI152_T1_1mm.nii.gz
ln -sf ../MNI152_T1_1mm_brain.nii.gz
ln -sf ../MNI152_T1_3mm.nii.gz

cd ../
for side in L R ; do
    for suffix in HEAD BRIK.gz ; do
	ln -sf sgACC/${side}.sgACC.mask.3mm+tlrc.${suffix}
    done
done

ln -sf sgACC/sgACC.mask.3mm+tlrc.HEAD
ln -sf sgACC/sgACC.mask.3mm+tlrc.BRIK.gz

