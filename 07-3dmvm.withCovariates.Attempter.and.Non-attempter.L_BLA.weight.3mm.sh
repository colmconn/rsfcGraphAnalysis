#!/bin/bash

### Variables in the data table are already mean centered

cd /data/sanDiego/rsfcGraphAnalysis/data//Group.results
3dMVM \
      -prefix restingstate.Attempter.and.Non-attempter.L_BLA.weight.3mm.3dmvm.bucket.20150825-1449PDT -mask mask.grey.Attempter.and.Non-attempter.union.masked+tlrc.HEAD \
      -jobs 8 \
      -ranEff '~1' \
      -SS_type 3 \
      -qVars 'Full' \
      -bsVars 'Group' \
      -wsVars 'Gender,age.in.years' \
      -num_glt  <n> \
      -gltLabel 1 '< add GLT here>' \
      -gltLabel 2 '< add GLT here>' \
      -dataTable @/data/sanDiego/rsfcGraphAnalysis/data//Group.data/dataTable.Attempter.and.Non-attempter.L_BLA.weight.3mm.txt
