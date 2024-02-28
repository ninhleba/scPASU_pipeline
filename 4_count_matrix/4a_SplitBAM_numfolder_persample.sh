#!/bin/bash

declare -A sample_numfolder

sample_numfolder=( ["BL_190994T"]="42" ["BL_191764T"]="74" ["BL_198490T"]="53" ["BL_199764"]="96" ["BL_207161"]="65" ["BL_208196"]="68" )  

#sample="BL_190994T BL_191764T BL_198490T BL_199764 BL_207161 BL_208196"
#sample="BL_190994T"
sample="BL_191764T BL_198490T BL_199764 BL_207161 BL_208196"

for s in $(echo $sample)
   do
     export n=${sample_numfolder[$s]}
     echo $s $n
     sed "s/samplename/$s/g" 4a_SplitBAM_parallel.lsf | sed "s/numberfolder/$n/g" | bsub
done
