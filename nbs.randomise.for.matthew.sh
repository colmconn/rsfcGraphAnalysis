#!/bin/bash

cd ../data/vbm.subject.list.from.matthew.n114/stats

for ii in 1 2 3 ; do
    randomise_parallel -i GM_mod_merg_s2 -m network${ii}.mask -o network${ii} -d design.mat -t design.con -T -n 10000
done
